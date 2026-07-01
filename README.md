# rsomics-vcf-indv-stats

Per-individual VCF statistics: TsTv substitution-type summary, singleton/doubleton allele sites, and per-sample mean depth — reimplementing vcftools `--TsTv-summary`, `--singletons`, and `--depth`.

## Usage

```
rsomics-vcf-indv-stats --mode <MODE> <VCF>

Modes:
  tstv-summary   Substitution-type counts + Ts/Tv totals (.TsTv.summary)
  singletons     Sites with 1 or 2 ALT allele copies (.singletons)
  depth          Per-sample mean FORMAT/DP (.idepth)
```

## Install

```
cargo install rsomics-vcf-indv-stats
```

## Origin

This crate is an independent Rust reimplementation of vcftools 0.1.17 based on:
- The vcftools documentation and man page
- Public file-format specifications (VCFv4.1/4.2)
- Black-box behaviour testing against the vcftools 0.1.17 binary

No source code from the LGPL vcftools upstream was used as reference during
implementation. Test fixtures are independently generated.

License: MIT OR Apache-2.0  
Upstream credit: vcftools <https://vcftools.github.io> (LGPL-3.0)
