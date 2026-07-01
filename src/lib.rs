//! Per-individual VCF statistics (vcftools --TsTv-summary / --singletons /
//! --depth).
//!
//! Three independent accumulators share one streaming VCF pass per mode:
//!
//! * `TsTvSummary` — biallelic-SNP substitution-type counts + Ts/Tv totals.
//! * `Singletons`  — sites where an ALT allele appears 1 or 2 times total.
//! * `Depth`       — per-sample mean FORMAT/DP.

pub mod vcf;

use std::path::Path;

use rsomics_common::Result;
use serde::Serialize;

// ── TsTv summary ────────────────────────────────────────────────────────────

/// Counts for each substitution model, in vcftools output order.
#[derive(Debug, Default, Clone, Serialize)]
pub struct TsTvSummary {
    pub ac: u64,
    pub ag: u64,
    pub at: u64,
    pub cg: u64,
    pub ct: u64,
    pub gt: u64,
}

impl TsTvSummary {
    pub fn ts(&self) -> u64 {
        self.ag + self.ct
    }

    pub fn tv(&self) -> u64 {
        self.ac + self.at + self.cg + self.gt
    }

    /// Render as the `.TsTv.summary` table vcftools emits.
    #[must_use]
    pub fn to_text(&self) -> String {
        format!(
            "MODEL\tCOUNT\nAC\t{}\nAG\t{}\nAT\t{}\nCG\t{}\nCT\t{}\nGT\t{}\nTs\t{}\nTv\t{}\n",
            self.ac,
            self.ag,
            self.at,
            self.cg,
            self.ct,
            self.gt,
            self.ts(),
            self.tv(),
        )
    }
}

/// Classify one biallelic SNP pair → update the corresponding counter.
///
/// `ref_base` and `alt_base` are ASCII uppercase bytes in `ACGT`. Returns
/// without updating when either byte is not in that alphabet.
pub fn classify_snp(stats: &mut TsTvSummary, ref_base: u8, alt_base: u8) {
    let (lo, hi) = if ref_base < alt_base {
        (ref_base, alt_base)
    } else {
        (alt_base, ref_base)
    };
    match (lo, hi) {
        (b'A', b'C') => stats.ac += 1,
        (b'A', b'G') => stats.ag += 1,
        (b'A', b'T') => stats.at += 1,
        (b'C', b'G') => stats.cg += 1,
        (b'C', b'T') => stats.ct += 1,
        (b'G', b'T') => stats.gt += 1,
        _ => {}
    }
}

// ── Singletons ──────────────────────────────────────────────────────────────

/// One singleton or doubleton site.
#[derive(Debug, Clone, Serialize)]
pub struct SingletonRow {
    pub chrom: String,
    pub pos: u64,
    pub kind: SingletonKind,
    pub allele: String,
    pub indv: String,
}

/// Whether a site is a singleton (S, 1 copy) or doubleton (D, 2 copies).
#[derive(Debug, Clone, Copy, Serialize)]
pub enum SingletonKind {
    S,
    D,
}

impl SingletonKind {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::S => "S",
            Self::D => "D",
        }
    }
}

/// The complete singletons table.
#[derive(Debug, Default, Clone, Serialize)]
pub struct Singletons {
    pub rows: Vec<SingletonRow>,
}

impl Singletons {
    /// Render as the `.singletons` table vcftools emits.
    #[must_use]
    pub fn to_text(&self) -> String {
        let mut out = String::from("CHROM\tPOS\tSINGLETON/DOUBLETON\tALLELE\tINDV\n");
        for r in &self.rows {
            out.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\n",
                r.chrom,
                r.pos,
                r.kind.as_str(),
                r.allele,
                r.indv,
            ));
        }
        out
    }
}

// ── Depth ────────────────────────────────────────────────────────────────────

/// Per-sample depth accumulator.
#[derive(Debug, Clone, Serialize)]
pub struct SampleDepth {
    pub sample: String,
    pub n_sites: u64,
    pub sum_dp: f64,
}

impl SampleDepth {
    fn mean(&self) -> f64 {
        if self.n_sites == 0 {
            0.0
        } else {
            self.sum_dp / self.n_sites as f64
        }
    }
}

/// Per-individual depth table.
#[derive(Debug, Default, Clone, Serialize)]
pub struct DepthTable {
    pub samples: Vec<SampleDepth>,
}

impl DepthTable {
    /// Render as the `.idepth` table vcftools emits.
    #[must_use]
    pub fn to_text(&self) -> String {
        let mut out = String::from("INDV\tN_SITES\tMEAN_DEPTH\n");
        for s in &self.samples {
            out.push_str(&format!(
                "{}\t{}\t{}\n",
                s.sample,
                s.n_sites,
                format_g6(s.mean()),
            ));
        }
        out
    }
}

/// Format a float with 6 significant figures, matching C `%g` (6 sig-figs, no
/// trailing zeros, no trailing decimal point).
pub fn format_g6(x: f64) -> String {
    if x == 0.0 {
        return "0".to_string();
    }
    // Determine how many integer digits the value has so we can compute the
    // number of decimal places needed for 6 total significant figures.
    let mag = x.abs().log10().floor() as i32;
    let dec = (5 - mag).max(0) as usize;
    let s = format!("{x:.dec$}");
    let s = s.trim_end_matches('0').trim_end_matches('.');
    s.to_string()
}

// ── Public entry points ──────────────────────────────────────────────────────

pub fn run_tstv_summary(path: &Path) -> Result<TsTvSummary> {
    vcf::scan_tstv(path)
}

pub fn run_singletons(path: &Path) -> Result<Singletons> {
    vcf::scan_singletons(path)
}

pub fn run_depth(path: &Path) -> Result<DepthTable> {
    vcf::scan_depth(path)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn classify_snp_all_models() {
        let mut s = TsTvSummary::default();
        classify_snp(&mut s, b'A', b'T'); // AT
        classify_snp(&mut s, b'G', b'C'); // CG (sorted)
        classify_snp(&mut s, b'C', b'T'); // CT
        assert_eq!(s.at, 1);
        assert_eq!(s.cg, 1);
        assert_eq!(s.ct, 1);
        assert_eq!(s.ts(), 1);
        assert_eq!(s.tv(), 2);
    }

    #[test]
    fn classify_snp_ignores_non_acgt() {
        let mut s = TsTvSummary::default();
        classify_snp(&mut s, b'N', b'A');
        assert_eq!(s.ts() + s.tv(), 0);
    }

    #[test]
    fn format_g6_integer() {
        assert_eq!(format_g6(18.0), "18");
    }

    #[test]
    fn format_g6_fractional() {
        // 32/3 = 10.666...
        let v = 32.0_f64 / 3.0;
        assert_eq!(format_g6(v), "10.6667");
    }

    #[test]
    fn format_g6_zero() {
        assert_eq!(format_g6(0.0), "0");
    }

    #[test]
    fn tstv_to_text() {
        let s = TsTvSummary {
            at: 1,
            cg: 1,
            ct: 1,
            ..Default::default()
        };
        let text = s.to_text();
        assert!(text.starts_with("MODEL\tCOUNT\n"));
        assert!(text.contains("AT\t1"));
        assert!(text.contains("Ts\t1"));
        assert!(text.contains("Tv\t2"));
    }

    #[test]
    fn singletons_to_text_empty() {
        let s = Singletons::default();
        assert_eq!(
            s.to_text(),
            "CHROM\tPOS\tSINGLETON/DOUBLETON\tALLELE\tINDV\n"
        );
    }

    #[test]
    fn depth_to_text() {
        let t = DepthTable {
            samples: vec![SampleDepth {
                sample: "S1".to_string(),
                n_sites: 3,
                sum_dp: 32.0,
            }],
        };
        let text = t.to_text();
        assert!(text.starts_with("INDV\tN_SITES\tMEAN_DEPTH\n"));
        assert!(text.contains("S1\t3\t10.6667"));
    }
}
