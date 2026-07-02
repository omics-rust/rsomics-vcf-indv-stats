# rsomics-vcf-indv-stats

Per-individual VCF statistics: TsTv substitution-type summary, singleton/doubleton allele sites, and per-sample mean depth — reimplementing vcftools `--TsTv-summary`, `--singletons`, and `--depth`.

## Usage

```
rsomics-vcf-indv-stats --mode <MODE> <VCF>

Modes:
  tstv-summary   Substitution-type counts + Ts/Tv totals (.TsTv.summary)
  singletons     Private singletons (S) and private doubletons (D) per allele (.singletons)
  depth          Per-sample mean FORMAT/DP (.idepth)
```

## Install

```
cargo install rsomics-vcf-indv-stats
```

## Boundaries

Output is byte-identical to vcftools 0.1.17 for every well-formed, coordinate-sorted,
equal-width VCF — including sites-only inputs, `POS=0`, missing/absent GT, half-calls,
phased and haploid genotypes, and multiallelic sites — and matches vcftools' exit code
and core stderr message on the common error cases (sites-only "Require Genotypes",
polyploid abort).

For genuinely malformed input, vcftools has undefined or buggy behaviour, so this crate
does **not** attempt to reproduce it and is instead the deterministic reference:

- **Ragged rows** — a data row with more (or fewer) sample columns than the `#CHROM`
  line declares. The header-declared sample count is authoritative; extra columns are
  ignored rather than indexed into a phantom sample.
- **Malformed GT separators** — tokens such as `.//` that split into more than two
  allele slots are treated as ploidy > 2 and abort like a polyploid site.
- **Duplicate FORMAT keys** — a FORMAT such as `GT:GT` resolves to the first `GT`
  slot; vcftools reads a different slot, so results differ.

## Origin

This crate is an independent Rust reimplementation of vcftools 0.1.17 based on:
- The vcftools documentation and man page
- Public file-format specifications (VCFv4.1/4.2)
- Black-box behaviour testing against the vcftools 0.1.17 binary

No source code from the LGPL vcftools upstream was used as reference during
implementation. Test fixtures are independently generated.

License: MIT OR Apache-2.0  
Upstream credit: vcftools <https://vcftools.github.io> (LGPL-3.0)
