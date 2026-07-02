//! Streaming byte-level VCF reader for per-individual statistics.
//!
//! Parses the text VCF format directly without a schema library so the hot
//! loop does no allocation per record beyond the initial sample-count sizing.
//! Supports plain and gzip-compressed inputs.

use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

use flate2::read::MultiGzDecoder;
use rsomics_common::{Result, RsomicsError};

use crate::{
    DepthTable, SampleDepth, SingletonKind, SingletonRow, SingletonScan, Singletons, TsTvSummary,
    classify_snp,
};

// ── I/O helper ───────────────────────────────────────────────────────────────

fn open_reader(path: &Path) -> Result<Box<dyn Read>> {
    let file = File::open(path).map_err(|e| {
        RsomicsError::Io(std::io::Error::new(
            e.kind(),
            format!("cannot open {}: {e}", path.display()),
        ))
    })?;
    let is_gz = path
        .extension()
        .is_some_and(|e| e.eq_ignore_ascii_case("gz"));
    let r: Box<dyn Read> = if is_gz {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };
    Ok(r)
}

fn missing_chrom_err() -> rsomics_common::RsomicsError {
    RsomicsError::InvalidInput("VCF missing #CHROM header line".into())
}

/// vcftools refuses `--singletons` on a genotype-less VCF (0 individuals),
/// exiting 1 with this exact core message.
fn require_genotypes_err() -> rsomics_common::RsomicsError {
    RsomicsError::InvalidInput(
        "Require Genotypes in VCF file in order to output Singletons.".into(),
    )
}

// ── Column indices (0-based after splitting on TAB) ──────────────────────────

const COL_CHROM: usize = 0;
const COL_POS: usize = 1;
const COL_REF: usize = 3;
const COL_ALT: usize = 4;
const COL_FORMAT: usize = 8;
const FIRST_SAMPLE: usize = 9;

/// Split a tab-delimited line into owned column strings. The trailing newline
/// has already been stripped by `BufRead::read_line`.
fn split_cols(line: &str) -> Vec<&str> {
    line.split('\t').collect()
}

// ── ACGT gate ────────────────────────────────────────────────────────────────

fn is_acgt(b: u8) -> bool {
    matches!(b, b'A' | b'C' | b'G' | b'T')
}

// ── TsTv summary ─────────────────────────────────────────────────────────────

/// Scan one VCF and accumulate biallelic-SNP substitution counts.
pub fn scan_tstv(path: &Path) -> Result<TsTvSummary> {
    let reader = open_reader(path)?;
    let mut lines = BufReader::new(reader).lines();
    let mut stats = TsTvSummary::default();
    let mut found_chrom = false;

    for line in lines.by_ref() {
        let line = line?;
        let line = line.trim_end_matches('\r');
        if line.starts_with("##") {
            continue;
        }
        if line.starts_with('#') {
            found_chrom = true;
            continue;
        }
        if !found_chrom {
            return Err(missing_chrom_err());
        }
        let cols = split_cols(line);
        if cols.len() < COL_ALT + 1 {
            continue;
        }
        let ref_col = cols[COL_REF].as_bytes();
        let alt_col = cols[COL_ALT];
        // Biallelic SNP only: REF is a single ACGT base, ALT has no comma.
        if ref_col.len() != 1 || !is_acgt(ref_col[0]) {
            continue;
        }
        if alt_col.contains(',') {
            continue;
        }
        let alt_bytes = alt_col.as_bytes();
        if alt_bytes.len() != 1 || !is_acgt(alt_bytes[0]) {
            continue;
        }
        classify_snp(&mut stats, ref_col[0], alt_bytes[0]);
    }

    if !found_chrom {
        return Err(missing_chrom_err());
    }
    Ok(stats)
}

// ── Singletons ───────────────────────────────────────────────────────────────

/// A diploid genotype call as vcftools stores it: two allele slots, each an
/// allele index or `-1` for a missing/absent copy. Haploid `1` → `(1, -1)`,
/// half-call `0/.` → `(0, -1)`, missing `./.` → `(-1, -1)`.
type GtPair = (i64, i64);

/// Parse a sample column's GT subfield into a `GtPair`. `gt_idx` is the GT
/// slot's position within the colon-delimited FORMAT; a sample lacking that
/// many subfields is treated as missing.
///
/// `Err` is reserved for ploidy > 2, which vcftools rejects outright
/// ("Polyploidy found, and not supported by vcftools").
fn parse_gt_pair(field: &str, gt_idx: usize) -> std::result::Result<GtPair, ()> {
    let gt = field.split(':').nth(gt_idx).unwrap_or(".");
    let mut alleles = gt.split(['/', '|']);
    let first = alleles.next();
    let second = alleles.next();
    if alleles.next().is_some() {
        return Err(());
    }
    let parse = |tok: Option<&str>| -> i64 {
        match tok {
            Some(t) if t != "." => t.parse::<i64>().unwrap_or(-1),
            _ => -1,
        }
    };
    Ok((parse(first), parse(second)))
}

/// Scan one VCF and collect singleton/doubleton rows.
///
/// For every allele index at a site — `0` = REF, `i` = the *i*-th ALT — whose
/// total non-missing copy count across all samples is exactly one or two,
/// vcftools emits a row: `S` for a single carrier, `D` only when both copies of
/// a doubleton sit in one homozygous individual. A doubleton split across two
/// carriers produces no row.
pub fn scan_singletons(path: &Path) -> Result<SingletonScan> {
    let reader = open_reader(path)?;
    let mut lines_iter = BufReader::new(reader).lines();
    let mut singletons = Singletons::default();
    let mut sample_names: Vec<String> = Vec::new();
    let mut found_chrom = false;
    let mut abort: Option<String> = None;

    for line in lines_iter.by_ref() {
        let line = line?;
        let line = line.trim_end_matches('\r');
        if line.starts_with("##") {
            continue;
        }
        if line.starts_with('#') {
            found_chrom = true;
            let cols: Vec<&str> = line.split('\t').collect();
            if cols.len() <= FIRST_SAMPLE {
                return Err(require_genotypes_err());
            }
            sample_names = cols[FIRST_SAMPLE..].iter().map(|s| s.to_string()).collect();
            continue;
        }
        if !found_chrom {
            return Err(missing_chrom_err());
        }
        let cols = split_cols(line);
        if cols.len() <= COL_FORMAT {
            continue;
        }
        // The GT slot's position inside FORMAT is per-site; a FORMAT without a
        // GT key carries no genotype data, so the whole site emits no rows.
        let gt_idx = match cols[COL_FORMAT].split(':').position(|k| k == "GT") {
            Some(i) => i,
            None => continue,
        };
        let chrom = cols[COL_CHROM];
        let pos: u64 = cols[COL_POS].parse().unwrap_or(0);
        let ref_col = cols[COL_REF];
        let alt_col = cols[COL_ALT];

        // Allele strings in vcftools index order: 0 = REF, then each ALT.
        let mut alleles: Vec<&str> = Vec::with_capacity(1 + alt_col.len());
        alleles.push(ref_col);
        if alt_col != "." {
            alleles.extend(alt_col.split(','));
        }
        let n_alleles = alleles.len();

        // vcftools reads every genotype for a site before emitting rows, so a
        // polyploid call aborts the whole site (no row for it) and terminates
        // the scan; rows from earlier sites have already been written.
        // The #CHROM header declares the sample count; a data row that carries
        // extra columns (ragged input) is truncated to that width so a sample
        // index can never outrun `sample_names`.
        let mut gts: Vec<GtPair> = Vec::with_capacity(sample_names.len());
        for field in cols[FIRST_SAMPLE..].iter().take(sample_names.len()) {
            match parse_gt_pair(field, gt_idx) {
                Ok(p) => gts.push(p),
                Err(()) => {
                    abort = Some(format!("{chrom}:{pos}"));
                    break;
                }
            }
        }
        if abort.is_some() {
            break;
        }

        let mut counts = vec![0u32; n_alleles];
        for &(a0, a1) in &gts {
            for a in [a0, a1] {
                if a >= 0 && (a as usize) < n_alleles {
                    counts[a as usize] += 1;
                }
            }
        }

        for ui in 0..n_alleles {
            let (kind, indv) = match counts[ui] {
                1 => {
                    let idx = gts
                        .iter()
                        .position(|&(a0, a1)| a0 == ui as i64 || a1 == ui as i64);
                    (SingletonKind::S, idx)
                }
                2 => {
                    let idx = gts
                        .iter()
                        .position(|&(a0, a1)| a0 == ui as i64 && a1 == ui as i64);
                    (SingletonKind::D, idx)
                }
                _ => continue,
            };
            let Some(sample_idx) = indv else { continue };
            singletons.rows.push(SingletonRow {
                chrom: chrom.to_string(),
                pos,
                kind,
                allele: alleles[ui].to_ascii_uppercase(),
                indv: sample_names[sample_idx].clone(),
            });
        }
    }

    if !found_chrom {
        return Err(missing_chrom_err());
    }
    Ok(SingletonScan { singletons, abort })
}

// ── Depth ─────────────────────────────────────────────────────────────────────

/// Scan one VCF and compute per-sample mean FORMAT/DP.
pub fn scan_depth(path: &Path) -> Result<DepthTable> {
    let reader = open_reader(path)?;
    let mut lines_iter = BufReader::new(reader).lines();
    let mut table = DepthTable::default();
    let mut found_chrom = false;

    for line in lines_iter.by_ref() {
        let line = line?;
        let line = line.trim_end_matches('\r');
        if line.starts_with("##") {
            continue;
        }
        if line.starts_with('#') {
            found_chrom = true;
            let cols: Vec<&str> = line.split('\t').collect();
            table.samples = cols[FIRST_SAMPLE..]
                .iter()
                .map(|n| SampleDepth {
                    sample: n.to_string(),
                    n_sites: 0,
                    sum_dp: 0.0,
                })
                .collect();
            continue;
        }
        if !found_chrom {
            return Err(missing_chrom_err());
        }
        if table.samples.is_empty() {
            continue;
        }
        let cols = split_cols(line);
        if cols.len() <= COL_FORMAT {
            continue;
        }
        // Find DP key index in FORMAT.
        let format_keys: Vec<&str> = cols[COL_FORMAT].split(':').collect();
        let dp_pos = match format_keys.iter().position(|&k| k == "DP") {
            Some(p) => p,
            None => continue,
        };

        for (sample_acc, sample_col) in table.samples.iter_mut().zip(cols[FIRST_SAMPLE..].iter()) {
            let fields: Vec<&str> = sample_col.split(':').collect();
            let dp_str = match fields.get(dp_pos) {
                Some(&s) => s,
                None => continue,
            };
            if dp_str == "." {
                continue;
            }
            let dp: f64 = match dp_str.parse() {
                Ok(v) => v,
                Err(_) => continue,
            };
            sample_acc.n_sites += 1;
            sample_acc.sum_dp += dp;
        }
    }

    if !found_chrom {
        return Err(missing_chrom_err());
    }
    Ok(table)
}
