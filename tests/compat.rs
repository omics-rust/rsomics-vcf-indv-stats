//! Value-exact compatibility with vcftools 0.1.17 per-individual statistics.
//!
//! All golden VCF content is hardcoded inline; no filesystem paths or external
//! processes are required. A second section gates on vcftools being on PATH and
//! diffs the live oracle output byte-for-byte.

use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

/// Three-sample VCF used for all three mode tests.
const GOLDEN_VCF: &str = "\
##fileformat=VCFv4.1\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\tSample2\tSample3\n\
chr1\t100\t.\tA\tT\t50\tPASS\t.\tGT:DP\t0/0:10\t0/1:15\t1/1:20\n\
chr1\t200\t.\tG\tC\t60\tPASS\t.\tGT:DP\t0/1:12\t1/1:18\t0/1:14\n\
chr2\t100\t.\tC\tT\t55\tPASS\t.\tGT:DP\t0/0:10\t0/0:10\t0/1:20\n\
";

/// Write the inline VCF bytes to a temp file in the KIOXIA tmp dir and return
/// the path. Uses `std::env::temp_dir()` so CI uses the runner's TMPDIR.
fn write_vcf(vcf: &str) -> PathBuf {
    let dir = std::env::temp_dir();
    let path = dir.join("rsomics_vcf_indv_stats_golden.vcf");
    let mut f = std::fs::File::create(&path).expect("cannot create temp VCF");
    f.write_all(vcf.as_bytes()).expect("write");
    path
}

// ── Expected outputs frozen from vcftools 0.1.17 black-box ───────────────────

const EXPECTED_TSTV: &str = "\
MODEL\tCOUNT\n\
AC\t0\n\
AG\t0\n\
AT\t1\n\
CG\t1\n\
CT\t1\n\
GT\t0\n\
Ts\t1\n\
Tv\t2\n\
";

const EXPECTED_SINGLETONS: &str = "\
CHROM\tPOS\tSINGLETON/DOUBLETON\tALLELE\tINDV\n\
chr2\t100\tS\tT\tSample3\n\
";

const EXPECTED_DEPTH: &str = "\
INDV\tN_SITES\tMEAN_DEPTH\n\
Sample1\t3\t10.6667\n\
Sample2\t3\t14.3333\n\
Sample3\t3\t18\n\
";

// ── Unit assertions against hardcoded expected strings ────────────────────────

#[test]
fn tstv_summary_matches_expected() {
    let path = write_vcf(GOLDEN_VCF);
    let stats = rsomics_vcf_indv_stats::run_tstv_summary(&path).unwrap();
    assert_eq!(stats.to_text(), EXPECTED_TSTV);
}

#[test]
fn singletons_matches_expected() {
    let path = write_vcf(GOLDEN_VCF);
    let singletons = rsomics_vcf_indv_stats::run_singletons(&path).unwrap();
    assert_eq!(singletons.to_text(), EXPECTED_SINGLETONS);
}

#[test]
fn depth_matches_expected() {
    let path = write_vcf(GOLDEN_VCF);
    let table = rsomics_vcf_indv_stats::run_depth(&path).unwrap();
    assert_eq!(table.to_text(), EXPECTED_DEPTH);
}

// ── Live vcftools oracle (skipped when vcftools not installed or wrong version)

fn vcftools_version() -> Option<String> {
    let out = Command::new("vcftools").arg("--version").output().ok()?;
    // vcftools prints version to stdout or stderr depending on build; check both.
    let combined =
        String::from_utf8_lossy(&out.stdout).to_string() + &String::from_utf8_lossy(&out.stderr);
    combined.lines().next().map(str::to_string)
}

fn skip_unless_vcftools_017() -> Option<()> {
    let ver = vcftools_version()?;
    if !ver.contains("0.1.17") {
        return None;
    }
    Some(())
}

/// Run vcftools with `--TsTv-summary` on a temp VCF and return the table text.
fn oracle_tstv(vcf: &Path) -> Option<String> {
    let out_dir = std::env::temp_dir();
    let prefix = out_dir.join("rsomics_vcf_indv_stats_oracle");
    let status = Command::new("vcftools")
        .args([
            "--vcf",
            vcf.to_str()?,
            "--TsTv-summary",
            "--out",
            prefix.to_str()?,
        ])
        .status()
        .ok()?;
    if !status.success() {
        return None;
    }
    let summary = prefix.with_extension("TsTv.summary");
    std::fs::read_to_string(summary).ok()
}

fn oracle_singletons(vcf: &Path) -> Option<String> {
    let out_dir = std::env::temp_dir();
    let prefix = out_dir.join("rsomics_vcf_indv_stats_oracle");
    let status = Command::new("vcftools")
        .args([
            "--vcf",
            vcf.to_str()?,
            "--singletons",
            "--out",
            prefix.to_str()?,
        ])
        .status()
        .ok()?;
    if !status.success() {
        return None;
    }
    let singletons = prefix.with_extension("singletons");
    std::fs::read_to_string(singletons).ok()
}

fn oracle_depth(vcf: &Path) -> Option<String> {
    let out_dir = std::env::temp_dir();
    let prefix = out_dir.join("rsomics_vcf_indv_stats_oracle");
    let status = Command::new("vcftools")
        .args(["--vcf", vcf.to_str()?, "--depth", "--out", prefix.to_str()?])
        .status()
        .ok()?;
    if !status.success() {
        return None;
    }
    let depth = prefix.with_extension("idepth");
    std::fs::read_to_string(depth).ok()
}

#[test]
fn live_oracle_tstv_summary() {
    if skip_unless_vcftools_017().is_none() {
        eprintln!("vcftools 0.1.17 not found — skipping live oracle TsTv test");
        return;
    }
    let path = write_vcf(GOLDEN_VCF);
    let oracle = oracle_tstv(&path).expect("vcftools --TsTv-summary failed");
    let stats = rsomics_vcf_indv_stats::run_tstv_summary(&path).unwrap();
    assert_eq!(stats.to_text(), oracle, "TsTv-summary differs from oracle");
}

#[test]
fn live_oracle_singletons() {
    if skip_unless_vcftools_017().is_none() {
        eprintln!("vcftools 0.1.17 not found — skipping live oracle singletons test");
        return;
    }
    let path = write_vcf(GOLDEN_VCF);
    let oracle = oracle_singletons(&path).expect("vcftools --singletons failed");
    let singletons = rsomics_vcf_indv_stats::run_singletons(&path).unwrap();
    assert_eq!(
        singletons.to_text(),
        oracle,
        "singletons differs from oracle"
    );
}

#[test]
fn live_oracle_depth() {
    if skip_unless_vcftools_017().is_none() {
        eprintln!("vcftools 0.1.17 not found — skipping live oracle depth test");
        return;
    }
    let path = write_vcf(GOLDEN_VCF);
    let oracle = oracle_depth(&path).expect("vcftools --depth failed");
    let table = rsomics_vcf_indv_stats::run_depth(&path).unwrap();
    assert_eq!(table.to_text(), oracle, "depth differs from oracle");
}
