# ANCESTRY_REF — 1000 Genomes-anchored ancestry PCs

Produces the ancestry principal components used as confounder covariates by
every QTL pipeline in this repo. Principal components are computed in the
1000 Genomes phase 3 reference **only**, then each cohort is projected into
that fixed space.

Everything runs locally — no GCP, no cluster.

```
ancestry/          the pipeline (Snakemake)
  Snakefile        13 rules, fetch -> merge -> shared SNP set -> PCA -> project -> report
  config.yaml      all paths, thresholds and cohort wiring; start here
  run.sh           wrapper: conda activate pd-qtl, snakemake --sdm conda
  scripts/         ancestry_projection.R (labels, plots, QC), cohort_keep_list.sh
  resources/       long_range_ld_hg19.bed — regions excluded before pruning
```

Build is GRCh37/hg19 with bare contig naming (`1`, not `chr1`) throughout.

## Running it

```bash
cd ancestry
./run.sh                 # CORES=4 by default; CORES=16 ./run.sh to go wider
./run.sh -n              # dry run
```

Safe to re-run: Snakemake resumes where it stopped, and the per-chromosome
1000G VCFs are `temp()` so they are deleted after merging rather than
re-downloaded.

## Outputs

Written to `results/ANCESTRY_REF/`:

| Path | Contents |
|---|---|
| `covariates/ancestry_pcs_PPMI.csv` | 894 samples — `IID,PC1..PC10` |
| `covariates/ancestry_pcs_PEG.csv` | 760 samples — `IID,PC1..PC10` |
| `ancestry_pcs_all.csv` | both cohorts plus the projected reference |
| `plots/ancestry_pca_cohorts.png` | cohorts over the 1000G background |
| `plots/ancestry_pca_facets.png` | per-superpopulation facets |
| `plots/ancestry_scree.png` | variance explained — the elbow is at PC5 |
| `qc/snp_set_summary.txt` | SNP counts surviving each filter |
| `qc/self_report_vs_pca.txt` | PCA result vs self-reported race/ethnicity |
| `qc/projection_check.txt` | reference-projected-through-its-own-weights check |
| `qc/peg_sample_selection.txt` | how 1840 VCF samples resolve to 760 PEG |

20 PCs are computed, 10 are written out, and the QTL pipelines currently use
**5** (`n_pcs: 5` in each of their configs) — matching the scree elbow.

Consumed by `PPMI_MQTL`, `PEG_MQTL`, `PEG_METAB` (via `config.yaml`) and
`PPMI_LONG` (via `run_longitudinal_mqtl.R`). These PCs are intended to
replace PEG's legacy STRUCTURE-based `K1a`–`K4a` admixture proportions.

## Why it is built this way

Four decisions here are deliberate and easy to undo by accident.

**Projection, not joint PCA.** PPMI is WGS and PEG is SNP array. A joint PCA
over the merged samples lets that platform difference — different
ascertainment, missingness, and error modes — load straight onto the top PCs,
which is precisely what a covariate meant to capture *ancestry* must not
contain. With projection the axes are fixed by 2504 reference samples, so
neither cohort can move them, and PPMI and PEG PCs end up in the same
coordinate system. That shared system is what makes the cross-cohort
meta-analysis in `META_MQTL` valid.

**The anchor is 1000G's released genotypes, not genotypes called from its
BAMs.** The phase 3 BAMs are ~4x low-coverage; the release is imputed and
refined across all 2504 samples and is a far better reference than anything
callable from the reads.

**Cohort inputs are pre-imputation genotypes.** PPMI's pre-imputation WGS
calls and PEG's pre-imputation array calls — deliberately not either cohort's
post-imputation dosages. Those dosages were themselves imputed against a
1000G-family panel, so projecting them onto a 1000G anchor partially
circularises the analysis and shrinks both cohorts artificially toward the
reference.

**Strand-ambiguous SNPs (A/T, C/G) are dropped.** Their strand cannot be
resolved by allele matching, and PEG is array data where a flip is entirely
plausible. A flipped ambiguous SNP becomes a silent allele-frequency error
that loads directly onto the top PCs.

Alongside those: variant IDs are rewritten to `CHROM:POS:REF:ALT` everywhere,
so an ID match also guarantees the alleles agree and REF/ALT-swapped variants
drop out instead of merging backwards; and long-range LD regions (the chr6
MHC and chr8p23 inversion) are excluded before pruning, since their internal
LD produces PCs reflecting local haplotype structure rather than genome-wide
ancestry.

## Traps worth knowing before you touch it

**PPMI input file.** Use `merged.dedup.vcf.gz`, never the
`plink2/all_chroms.*` fileset. That fileset predates the sample-duplication
fix and carries every patient twice — as `PPMISI{id}.variant` (an indel
callset that is a near-constant `0/0` stub at SNP sites) and
`PPMISI{id}.variant2` (the real SNP genotype). Running the PCA on it halves
every allele frequency and plots each patient twice.

**PEG sample identity.** The PEG VCF's sample names are GWAS_IDs
(`CRG_*`/`CRG2_*`), while the covariate file is keyed by `Pegid`.
Intersecting the two namespaces directly yields the empty set. They are
joined through `gwas_linkage.csv`, dropping `DUP=1` rows so the join is
exactly 1:1 — the same file and convention `PEG_MQTL` and `PEG_METAB` use.

**`snp_set.ld_prune.step` must be 1.** plink2 accepts a window increment > 1
only when the window is counted in variants, and this window is in kb. The
Snakefile asserts this at parse time rather than letting a run fail 20
minutes in.

**The shared SNP set is small.** 6.86 M common 1000G variants intersected
with PPMI (33.7 M) and PEG (304 K) leaves 38,329 before pruning — PEG's array
content is the binding constraint. Check `qc/snp_set_summary.txt` after any
change to the filters; a sharp drop there is the first sign something stopped
matching.

## Not part of this pipeline

The read-level NGS exercise that used to live here (`ngs_exercise/`) has moved
to its own repository. It was a self-contained learning exercise — aligning
and variant-calling GIAB reads to benchmark a workflow — and it never fed
these ancestry PCs or any QTL covariate. Keep it that way.
