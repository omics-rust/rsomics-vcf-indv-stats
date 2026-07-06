//! Value-exact compatibility with vcftools 0.1.17 per-individual statistics.
//!
//! Every expected string is frozen from a live `vcftools 0.1.17 --singletons /
//! --TsTv-summary / --depth` run captured once and pasted here as a constant.
//! The tests never spawn vcftools, Python, or any external process; they only
//! write a tiny VCF to the process temp dir and read it back through the crate.

use std::io::Write;
use std::path::PathBuf;

/// Three-sample VCF used for the TsTv/depth mode tests and the baseline
/// singletons case.
const GOLDEN_VCF: &str = "\
##fileformat=VCFv4.1\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\tSample2\tSample3\n\
chr1\t100\t.\tA\tT\t50\tPASS\t.\tGT:DP\t0/0:10\t0/1:15\t1/1:20\n\
chr1\t200\t.\tG\tC\t60\tPASS\t.\tGT:DP\t0/1:12\t1/1:18\t0/1:14\n\
chr2\t100\t.\tC\tT\t55\tPASS\t.\tGT:DP\t0/0:10\t0/0:10\t0/1:20\n\
";

/// Every singleton/doubleton edge in one file: spread doubleton (no row),
/// private doubleton, REF singleton, ALT singleton, half-call carrying ALT,
/// half-call carrying REF, phased, haploid spread doubleton (no row → REF
/// singleton), multiallelic double singleton, monomorphic, all-missing, indel,
/// symbolic ALT.
const EDGES_VCF: &str = "\
##fileformat=VCFv4.1\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n\
chr1\t10\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/1\t0/0\n\
chr1\t20\t.\tA\tT\t.\t.\t.\tGT\t1/1\t0/0\t0/0\n\
chr1\t30\t.\tA\tT\t.\t.\t.\tGT\t1/1\t1/1\t0/1\n\
chr1\t40\t.\tA\tT\t.\t.\t.\tGT\t0/0\t0/0\t0/1\n\
chr1\t50\t.\tA\tT\t.\t.\t.\tGT\t0/0\t0/0\t./1\n\
chr1\t60\t.\tA\tT\t.\t.\t.\tGT\t1/1\t1/1\t0/.\n\
chr1\t70\t.\tA\tT\t.\t.\t.\tGT\t0|0\t0|0\t0|1\n\
chr1\t80\t.\tA\tT\t.\t.\t.\tGT\t1\t1\t0\n\
chr1\t90\t.\tA\tT,G\t.\t.\t.\tGT\t0/0\t0/1\t0/2\n\
chr1\t100\t.\tA\tT\t.\t.\t.\tGT\t0/0\t0/0\t0/0\n\
chr1\t110\t.\tA\tT\t.\t.\t.\tGT\t./.\t./.\t./.\n\
chr2\t5\t.\tAC\tA\t.\t.\t.\tGT\t0/0\t0/0\t0/1\n\
chr2\t15\t.\tA\t<DEL>\t.\t.\t.\tGT\t0/0\t0/0\t0/1\n\
";

/// One sample, heterozygous 0/1: both REF and ALT have count one, so vcftools
/// emits a singleton row for each.
const SINGLE_SAMPLE_VCF: &str = "\
##fileformat=VCFv4.1\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tONLY\n\
chr1\t100\t.\tA\tT\t.\t.\t.\tGT\t0/1\n\
";

/// A polyploid genotype at a site — vcftools aborts with a hard error.
const POLYPLOID_VCF: &str = "\
##fileformat=VCFv4.1\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n\
chr1\t100\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/0/0\t0/0\n\
";

/// Write VCF bytes to a uniquely-named temp file and return the path. Each call
/// produces a distinct name to avoid races between parallel tests.
fn write_vcf(vcf: &str) -> PathBuf {
    use std::time::{SystemTime, UNIX_EPOCH};
    let dir = std::env::temp_dir();
    let ts = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.subsec_nanos())
        .unwrap_or(0);
    let tid = std::thread::current().id();
    let name = format!("rsomics_vcf_indv_stats_{tid:?}_{ts}.vcf");
    let path = dir.join(name);
    let mut f = std::fs::File::create(&path).expect("cannot create temp VCF");
    f.write_all(vcf.as_bytes()).expect("write");
    path
}

// ── Expected outputs frozen from vcftools 0.1.17 ─────────────────────────────

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

const EXPECTED_EDGES: &str = "\
CHROM\tPOS\tSINGLETON/DOUBLETON\tALLELE\tINDV\n\
chr1\t20\tD\tT\tS1\n\
chr1\t30\tS\tA\tS3\n\
chr1\t40\tS\tT\tS3\n\
chr1\t50\tS\tT\tS3\n\
chr1\t60\tS\tA\tS3\n\
chr1\t70\tS\tT\tS3\n\
chr1\t80\tS\tA\tS3\n\
chr1\t90\tS\tT\tS2\n\
chr1\t90\tS\tG\tS3\n\
chr2\t5\tS\tA\tS3\n\
chr2\t15\tS\t<DEL>\tS3\n\
";

const EXPECTED_SINGLE_SAMPLE: &str = "\
CHROM\tPOS\tSINGLETON/DOUBLETON\tALLELE\tINDV\n\
chr1\t100\tS\tA\tONLY\n\
chr1\t100\tS\tT\tONLY\n\
";

const EXPECTED_DEPTH: &str = "\
INDV\tN_SITES\tMEAN_DEPTH\n\
Sample1\t3\t10.6667\n\
Sample2\t3\t14.3333\n\
Sample3\t3\t18\n\
";

// ── Assertions against the frozen vcftools 0.1.17 output ─────────────────────

#[test]
fn tstv_summary_matches_expected() {
    let path = write_vcf(GOLDEN_VCF);
    let stats = rsomics_vcf_indv_stats::run_tstv_summary(&path).unwrap();
    assert_eq!(stats.to_text(), EXPECTED_TSTV);
}

#[test]
fn singletons_baseline_matches_expected() {
    let path = write_vcf(GOLDEN_VCF);
    let scan = rsomics_vcf_indv_stats::run_singletons(&path).unwrap();
    assert!(scan.abort.is_none());
    assert_eq!(scan.singletons.to_text(), EXPECTED_SINGLETONS);
}

#[test]
fn singletons_edges_match_expected() {
    let path = write_vcf(EDGES_VCF);
    let scan = rsomics_vcf_indv_stats::run_singletons(&path).unwrap();
    assert!(scan.abort.is_none());
    assert_eq!(scan.singletons.to_text(), EXPECTED_EDGES);
}

#[test]
fn singletons_single_sample_matches_expected() {
    let path = write_vcf(SINGLE_SAMPLE_VCF);
    let scan = rsomics_vcf_indv_stats::run_singletons(&path).unwrap();
    assert!(scan.abort.is_none());
    assert_eq!(scan.singletons.to_text(), EXPECTED_SINGLE_SAMPLE);
}

#[test]
fn singletons_small_frequency_matches_expected() {
    // 200 samples, one heterozygous ALT carrier → a single ALT singleton row.
    let mut vcf =
        String::from("##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    for i in 0..200 {
        vcf.push_str(&format!("\tN{i}"));
    }
    vcf.push_str("\nchr1\t100\t.\tA\tT\t.\t.\t.\tGT");
    for i in 0..200 {
        vcf.push_str(if i == 199 { "\t0/1" } else { "\t0/0" });
    }
    vcf.push('\n');
    let path = write_vcf(&vcf);
    let scan = rsomics_vcf_indv_stats::run_singletons(&path).unwrap();
    assert!(scan.abort.is_none());
    assert_eq!(
        scan.singletons.to_text(),
        "CHROM\tPOS\tSINGLETON/DOUBLETON\tALLELE\tINDV\nchr1\t100\tS\tT\tN199\n",
    );
}

#[test]
fn singletons_polyploid_is_rejected() {
    // vcftools aborts a --singletons run when a genotype has ploidy > 2, but
    // only after writing a header-only file (the polyploid is the first site).
    let path = write_vcf(POLYPLOID_VCF);
    let scan = rsomics_vcf_indv_stats::run_singletons(&path).unwrap();
    assert_eq!(scan.abort.as_deref(), Some("chr1:100"));
    assert_eq!(
        scan.singletons.to_text(),
        "CHROM\tPOS\tSINGLETON/DOUBLETON\tALLELE\tINDV\n",
    );
}

/// Lowercase REF/ALT and mixed-case symbolic alleles — vcftools uppercases the
/// emitted ALLELE column universally (`a`→`A`, `<del>`→`<DEL>`).
const LOWERCASE_VCF: &str = "\
##fileformat=VCFv4.1\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n\
1\t100\t.\tacgt\ta\t.\t.\t.\tGT\t0/0\t0/0\t0/1\n\
1\t200\t.\tn\ta\t.\t.\t.\tGT\t0/0\t0/0\t0/1\n\
1\t300\t.\tA\t<del>\t.\t.\t.\tGT\t0/0\t0/0\t0/1\n\
1\t500\t.\tc\tA\t.\t.\t.\tGT\t1/1\t0/0\t0/0\n\
";

const EXPECTED_LOWERCASE: &str = "\
CHROM\tPOS\tSINGLETON/DOUBLETON\tALLELE\tINDV\n\
1\t100\tS\tA\tS3\n\
1\t200\tS\tA\tS3\n\
1\t300\tS\t<DEL>\tS3\n\
1\t500\tD\tA\tS1\n\
";

#[test]
fn singletons_alleles_are_uppercased() {
    let path = write_vcf(LOWERCASE_VCF);
    let scan = rsomics_vcf_indv_stats::run_singletons(&path).unwrap();
    assert!(scan.abort.is_none());
    assert_eq!(scan.singletons.to_text(), EXPECTED_LOWERCASE);
}

/// A polyploid site preceded by real singleton rows: vcftools emits every row
/// before the offending site, then aborts. Frozen from vcftools 0.1.17.
const POLYPLOID_AFTER_ROWS_VCF: &str = "\
##fileformat=VCFv4.1\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n\
1\t100\t.\ta\tt\t.\t.\t.\tGT\t0/1\t0/0\n\
1\t200\t.\tA\tG\t.\t.\t.\tGT\t0/0/1\t0/0\n\
";

#[test]
fn singletons_polyploid_emits_preceding_rows() {
    let path = write_vcf(POLYPLOID_AFTER_ROWS_VCF);
    let scan = rsomics_vcf_indv_stats::run_singletons(&path).unwrap();
    assert_eq!(scan.abort.as_deref(), Some("1:200"));
    assert_eq!(
        scan.singletons.to_text(),
        "CHROM\tPOS\tSINGLETON/DOUBLETON\tALLELE\tINDV\n1\t100\tS\tT\tS1\n",
    );
}

/// FORMAT places GT after DP, so the genotype must be located by name, not by
/// taking the first colon-subfield. Frozen from vcftools 0.1.17.
const GT_NOT_FIRST_VCF: &str = "\
##fileformat=VCFv4.1\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n\
chr1\t100\t.\tA\tT\t.\t.\t.\tDP:GT\t10:0/0\t15:0/1\t20:0/0\n\
chr1\t200\t.\tG\tC\t.\t.\t.\tDP:GT\t12:0/0\t18:0/0\t14:0/1\n\
";

const EXPECTED_GT_NOT_FIRST: &str = "\
CHROM\tPOS\tSINGLETON/DOUBLETON\tALLELE\tINDV\n\
chr1\t100\tS\tT\tS2\n\
chr1\t200\tS\tC\tS3\n\
";

#[test]
fn singletons_gt_located_by_name_not_first_subfield() {
    let path = write_vcf(GT_NOT_FIRST_VCF);
    let scan = rsomics_vcf_indv_stats::run_singletons(&path).unwrap();
    assert!(scan.abort.is_none());
    assert_eq!(scan.singletons.to_text(), EXPECTED_GT_NOT_FIRST);
}

/// FORMAT carries no GT key though sample columns exist. vcftools keeps the 3
/// individuals but emits a header-only table — no site has genotype data.
const GT_ABSENT_VCF: &str = "\
##fileformat=VCFv4.1\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n\
chr1\t100\t.\tA\tT\t.\t.\t.\tDP\t1\t1\t1\n\
chr1\t200\t.\tG\tC\t.\t.\t.\tDP\t5\t0\t2\n\
";

#[test]
fn singletons_format_without_gt_emits_header_only() {
    let path = write_vcf(GT_ABSENT_VCF);
    let scan = rsomics_vcf_indv_stats::run_singletons(&path).unwrap();
    assert!(scan.abort.is_none());
    assert_eq!(
        scan.singletons.to_text(),
        "CHROM\tPOS\tSINGLETON/DOUBLETON\tALLELE\tINDV\n",
    );
}

/// A sites-only VCF (no genotype columns) makes vcftools exit 1 with a fixed
/// "Require Genotypes" message before reading any data.
const SITES_ONLY_VCF: &str = "\
##fileformat=VCFv4.1\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t100\t.\tA\tT\t.\t.\t.\n\
";

#[test]
fn singletons_sites_only_requires_genotypes() {
    let path = write_vcf(SITES_ONLY_VCF);
    let err = rsomics_vcf_indv_stats::run_singletons(&path).unwrap_err();
    assert!(
        err.to_string()
            .contains("Require Genotypes in VCF file in order to output Singletons."),
        "unexpected error: {err}"
    );
}

#[test]
fn depth_matches_expected() {
    let path = write_vcf(GOLDEN_VCF);
    let table = rsomics_vcf_indv_stats::run_depth(&path).unwrap();
    assert_eq!(table.to_text(), EXPECTED_DEPTH);
}

/// A sites-only VCF (8-column header, no FORMAT/samples) makes vcftools exit 1
/// with "...output Individuals by Mean Depth Statistics." — the depth analogue
/// of the singletons Require-Genotypes guard.
const DEPTH_SITES_ONLY_VCF: &str = "\
##fileformat=VCFv4.1\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t100\t.\tA\tT\t.\t.\t.\n\
";

/// A 9-column header (FORMAT present, zero sample columns) is likewise a
/// genotype-less file; vcftools exits 1 with the same message.
const DEPTH_ZERO_SAMPLE_VCF: &str = "\
##fileformat=VCFv4.1\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n\
chr1\t100\t.\tA\tT\t.\t.\t.\tGT:DP\n\
";

const DEPTH_REQUIRE_GT_MSG: &str =
    "Require Genotypes in VCF file in order to output Individuals by Mean Depth Statistics.";

#[test]
fn depth_sites_only_requires_genotypes() {
    let path = write_vcf(DEPTH_SITES_ONLY_VCF);
    let err = rsomics_vcf_indv_stats::run_depth(&path).unwrap_err();
    assert!(
        err.to_string().contains(DEPTH_REQUIRE_GT_MSG),
        "unexpected error: {err}"
    );
}

#[test]
fn depth_zero_sample_requires_genotypes() {
    let path = write_vcf(DEPTH_ZERO_SAMPLE_VCF);
    let err = rsomics_vcf_indv_stats::run_depth(&path).unwrap_err();
    assert!(
        err.to_string().contains(DEPTH_REQUIRE_GT_MSG),
        "unexpected error: {err}"
    );
}

/// The binary must exit non-zero and print the message to stderr, not panic,
/// on a genotype-less depth run.
#[test]
fn depth_fail_loud_exit_and_stderr() {
    for vcf in [DEPTH_SITES_ONLY_VCF, DEPTH_ZERO_SAMPLE_VCF] {
        let path = write_vcf(vcf);
        let out = std::process::Command::new(env!("CARGO_BIN_EXE_rsomics-vcf-indv-stats"))
            .args(["--mode", "depth"])
            .arg(&path)
            .output()
            .expect("spawn binary");
        assert!(!out.status.success(), "expected non-zero exit");
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert!(
            stderr.contains(DEPTH_REQUIRE_GT_MSG),
            "stderr missing require-genotypes message: {stderr}"
        );
        assert!(out.stdout.is_empty(), "expected no stdout on fail-loud");
    }
}

/// A sample whose FORMAT/DP is missing (`.`) at every site has N_SITES=0;
/// vcftools prints `nan` for the 0/0 mean. Frozen from vcftools 0.1.17.
const DEPTH_NAN_VCF: &str = "\
##fileformat=VCFv4.1\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n\
chr1\t100\t.\tA\tT\t50\tPASS\t.\tGT:DP\t0/0:10\t0/0:.\n\
chr1\t200\t.\tG\tC\t60\tPASS\t.\tGT:DP\t0/1:12\t0/1:.\n\
";

const EXPECTED_DEPTH_NAN: &str = "\
INDV\tN_SITES\tMEAN_DEPTH\n\
S1\t2\t11\n\
S2\t0\tnan\n\
";

#[test]
fn depth_zero_sites_prints_nan() {
    let path = write_vcf(DEPTH_NAN_VCF);
    let table = rsomics_vcf_indv_stats::run_depth(&path).unwrap();
    assert_eq!(table.to_text(), EXPECTED_DEPTH_NAN);
}

/// FORMAT/DP is spec-Integer; vcftools reads it with atoi, so `10.5` truncates
/// to `10` (mean 8.5, not 8.75). Frozen from vcftools 0.1.17.
const DEPTH_FLOAT_DP_VCF: &str = "\
##fileformat=VCFv4.1\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n\
chr1\t100\t.\tA\tT\t.\t.\t.\tGT:DP\t0/1:10.5\n\
chr1\t200\t.\tG\tC\t.\t.\t.\tGT:DP\t0/1:7\n\
";

const EXPECTED_DEPTH_FLOAT_DP: &str = "\
INDV\tN_SITES\tMEAN_DEPTH\n\
S1\t2\t8.5\n\
";

#[test]
fn depth_dp_truncated_to_integer() {
    let path = write_vcf(DEPTH_FLOAT_DP_VCF);
    let table = rsomics_vcf_indv_stats::run_depth(&path).unwrap();
    assert_eq!(table.to_text(), EXPECTED_DEPTH_FLOAT_DP);
}

/// Lowercase (soft-masked) biallelic-SNP bases are case-insensitive per the VCF
/// spec; vcftools folds them before the ACGT classification. Frozen from
/// vcftools 0.1.17: `a>g` is a transition, `c>t` is a transition → Ts=2.
const TSTV_LOWERCASE_VCF: &str = "\
##fileformat=VCFv4.1\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n\
chr1\t100\t.\ta\tg\t.\t.\t.\tGT\t0/1\n\
chr1\t200\t.\tc\tt\t.\t.\t.\tGT\t0/1\n\
";

const EXPECTED_TSTV_LOWERCASE: &str = "\
MODEL\tCOUNT\n\
AC\t0\n\
AG\t1\n\
AT\t0\n\
CG\t0\n\
CT\t1\n\
GT\t0\n\
Ts\t2\n\
Tv\t0\n\
";

#[test]
fn tstv_lowercase_bases_are_counted() {
    let path = write_vcf(TSTV_LOWERCASE_VCF);
    let stats = rsomics_vcf_indv_stats::run_tstv_summary(&path).unwrap();
    assert_eq!(stats.to_text(), EXPECTED_TSTV_LOWERCASE);
}
