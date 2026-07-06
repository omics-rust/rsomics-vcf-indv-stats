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
equal-width VCF — including `POS=0`, missing/absent GT, half-calls, phased and haploid
genotypes, multiallelic sites, lowercase/soft-masked REF/ALT bases (folded before the
ACGT gate, as vcftools does), non-integer FORMAT/DP (read with `atoi`: `10.5` → `10`),
and samples covered at zero sites (mean printed as `nan`).

It matches vcftools' exit code and core stderr message on the shared error cases:

- **Genotype-less `--singletons` or `--depth`** — a sites-only VCF (no sample columns),
  or a `#CHROM` line carrying a FORMAT column but zero samples, exits 1 with the
  respective "Require Genotypes in VCF file in order to output …" message.
- **Polyploid abort** — a genotype of ploidy > 2 in `--singletons` prints every row
  found before the offending site, then exits 1.

(`--TsTv-summary` needs no genotypes, so a sites-only VCF is valid input there and
produces a normal — possibly all-zero — table.)

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
