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
    DepthTable, SampleDepth, SingletonKind, SingletonRow, Singletons, TsTvSummary, classify_snp,
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

/// Count how many times each ALT allele index appears across all samples and
/// collect which samples carry it. Returns `(total_count, carriers)` per ALT.
///
/// `carriers[i]` = list of sample indices that carry at least one copy of
/// `ALT[i]` (0-based ALT index). `total[i]` = total allele copies.
fn count_alt_alleles(gt_fields: &[&str]) -> (Vec<u32>, Vec<Vec<usize>>) {
    // We don't know the ALT count ahead of time; size on demand.
    let mut totals: Vec<u32> = Vec::new();
    let mut carriers: Vec<Vec<usize>> = Vec::new();

    for (sample_idx, field) in gt_fields.iter().enumerate() {
        let gt = if let Some(colon) = field.find(':') {
            &field[..colon]
        } else {
            field
        };
        // Genotype alleles separated by '/' or '|'; skip missing '.'.
        let alleles: Vec<&str> = gt.split(['/', '|']).collect();
        for a in alleles {
            if a == "." {
                continue;
            }
            let idx: usize = match a.parse::<usize>() {
                Ok(v) if v > 0 => v - 1, // 0-based ALT index
                _ => continue,
            };
            // Extend vecs if needed.
            while totals.len() <= idx {
                totals.push(0);
                carriers.push(Vec::new());
            }
            totals[idx] += 1;
            let c = &mut carriers[idx];
            if c.last() != Some(&sample_idx) {
                c.push(sample_idx);
            }
        }
    }
    (totals, carriers)
}

/// Scan one VCF and collect singleton/doubleton rows.
pub fn scan_singletons(path: &Path) -> Result<Singletons> {
    let reader = open_reader(path)?;
    let mut lines_iter = BufReader::new(reader).lines();
    let mut singletons = Singletons::default();
    let mut sample_names: Vec<String> = Vec::new();
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
            sample_names = cols[FIRST_SAMPLE..].iter().map(|s| s.to_string()).collect();
            continue;
        }
        if !found_chrom {
            return Err(missing_chrom_err());
        }
        if sample_names.is_empty() {
            continue;
        }
        let cols = split_cols(line);
        if cols.len() < FIRST_SAMPLE {
            continue;
        }
        let chrom = cols[COL_CHROM].to_string();
        let pos: u64 = cols[COL_POS].parse().unwrap_or(0);
        let alt_col = cols[COL_ALT];
        let alts: Vec<&str> = alt_col.split(',').collect();
        let gt_fields = &cols[FIRST_SAMPLE..];

        let (totals, carriers) = count_alt_alleles(gt_fields);

        for (alt_idx, (&total, carrier_list)) in totals.iter().zip(carriers.iter()).enumerate() {
            let kind = match total {
                1 => SingletonKind::S,
                2 => SingletonKind::D,
                _ => continue,
            };
            let allele = alts.get(alt_idx).copied().unwrap_or(".").to_string();
            for &sample_idx in carrier_list {
                let indv = sample_names.get(sample_idx).cloned().unwrap_or_default();
                singletons.rows.push(SingletonRow {
                    chrom: chrom.clone(),
                    pos,
                    kind,
                    allele: allele.clone(),
                    indv,
                });
            }
        }
    }

    if !found_chrom {
        return Err(missing_chrom_err());
    }
    Ok(singletons)
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
