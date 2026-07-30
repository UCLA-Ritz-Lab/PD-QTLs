# PPMI WGS Pipeline: Per-Sample VCFs → Merged VCF → PLINK2

Processes per-sample whole-genome sequencing VCF files from the
[Parkinson's Progression Markers Initiative (PPMI)](https://www.ppmi-info.org/)
into a merged multi-sample VCF and PLINK2 binary files ready for mQTL analysis.

**Reference genome:** GRCh37 / hg19 (chromosomes named `1, 2, ... 22`)  
**Input:** Per-sample GATK-processed VCFs (~500+ samples)  
**Output:** PLINK2 `.pgen/.psam/.pvar` files, QC missingness report

---

## Pipeline Overview

```
Downloads/PPMISI{SAMPLE_ID}.realigned.recalibrated.g.genotyped.filtered.annotated.vcf.gz (500+)
    │
    ├─ Step 1: index_vcf          tabix index each input VCF (sentinel file)
    ├─ Step 2: normalize_vcf      strip malformed INFO tags, decompose
    │                             multiallelics, left-align indels, GQ/DP filter
    │
    └─ Step 3: merge_vcfs         bcftools merge → multi-sample VCF
                │
                ├─ Step 4: filter_merged     site-level missingness filter
                │
                └─ Step 5: split_by_chrom   one VCF per chromosome
                            │
                            └─ Step 6: vcf_to_plink2   per-chrom PLINK2 files
                                        │
                                        ├─ Step 7: merge_plink2   genome-wide pfiles
                                        └─ Step 8: qc_report      missingness summary
```

---

## Directory Structure

```
ppmi_wgs_pipeline/
├── Snakefile
├── config.yaml
├── README.md
├── envs/
│   ├── bcftools.yaml
│   ├── plink2.yaml
│   └── pyvcf.yaml
├── resources/
│   └── hg19.fa           ← download separately (see Setup)
│   └── hg19.fa.fai
└── Downloads/
    └── PPMISI{ID}.realigned.recalibrated.g.genotyped.filtered.annotated.vcf.gz
```

---

## Setup

### 1. Install Snakemake

```bash
conda create -n snakemake -c conda-forge -c bioconda snakemake
conda activate snakemake
```

Snakemake ≥ 7.0 is required. The `--use-conda` flag creates per-rule
environments automatically from the `envs/` YAML files.

### 2. Download the GRCh37 Reference Genome

Use the 1000 Genomes GRCh37 reference (`hs37d5`) — this is the correct
reference for PPMI WGS data as it includes decoy sequences (EBV/NC_007605
and others) that the PPMI GATK pipeline mapped to. Without these decoy
sequences, `bcftools norm` will fail on samples with reads mapped to
non-standard contigs.

```bash
mkdir -p resources

wget -P resources/ \
    ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz

bgzip -d resources/hs37d5.fa.gz
mv resources/hs37d5.fa resources/hg19.fa
samtools faidx resources/hg19.fa
```

> **Why not UCSC hg19?**  
> The UCSC hg19 reference uses `chr1` style chromosome names while PPMI
> VCFs use bare `1` style names. Using the wrong reference causes
> `bcftools norm --check-ref` to fail on every record. Always use the
> 1000 Genomes GRCh37 / hs37d5 reference for PPMI data.

### 3. Verify Chromosome Naming in Your VCFs

```bash
bcftools view -H \
    Downloads/PPMISI001.realigned.recalibrated.g.genotyped.filtered.annotated.vcf.gz \
    | head -3 | cut -f1
```

You should see `1`, `2`, `3` — not `chr1`. If you see `chr` prefixes,
update the `chromosomes` list in `config.yaml` accordingly and switch to
the UCSC hg19 reference.

### 4. Configure the Pipeline

Edit `config.yaml` to match your environment:

```yaml
vcf_dir:   "../../../Downloads/PPMI_WGS/pre_imputed"  # directory containing per-sample VCFs
out_dir:   "../../../results/PPMI_WGS/pre_imputed"    # all outputs go here
ref_fasta: "../../../Downloads/PPMI_WGS/resources/hg19.fa"  # path to hs37d5 reference
```

Adjust thread counts to match your workstation:

```yaml
resources:
  threads_bcftools: 4    # per normalisation job
  threads_merge: 8       # bcftools merge (most memory-intensive step)
  threads_plink: 8       # PLINK2 conversion and merge
  mem_gb_merge: 32       # raise if merge runs out of RAM (500+ samples)
```

### 5. Dry Run

Always do a dry run first to confirm sample discovery and rule counts:

```bash
conda activate snakemake
snakemake --use-conda --cores 1 --dry-run
```

The output should show one `normalize_vcf` job per sample. Confirm the
count matches your expected number of input VCF files.

### 6. Run the Pipeline

```bash
# Recommended: use all available cores
snakemake --use-conda --cores all

# Or set explicitly (e.g. 16-core workstation)
snakemake --use-conda --cores 16
```

For runs that will take many hours, use `tmux` or `screen` to prevent
disconnection from killing the process:

```bash
tmux new -s ppmi_wgs
snakemake --use-conda --cores all
# Detach: Ctrl+B then D
# Reattach: tmux attach -t ppmi_wgs
```

---

## Output Files

```
results/
├── index_sentinel/          completion sentinels for tabix indexing
├── normalised/              per-sample normalised VCFs (.norm.vcf.gz)
├── merged/
│   ├── vcf_list.txt                  list of input VCFs for bcftools merge
│   ├── merged.vcf.gz                 all samples merged, pre-filter
│   └── merged.filtered.vcf.gz        missingness-filtered  ← keep this
├── by_chrom/                per-chromosome VCFs (intermediate)
├── plink2/
│   ├── pmerge_list.txt               list of per-chrom pfiles for merge
│   ├── by_chrom/                     per-chromosome .pgen/.psam/.pvar
│   └── all_chroms.{pgen,psam,pvar}   ← final output for mQTL analysis
├── qc/
│   └── missingness_report.txt        per-sample and per-variant missingness
└── logs/                    per-rule log files for debugging
```

---


## PPMI-Specific Notes

### Malformed INFO Tags in PPMI VCFs

PPMI VCFs are annotated with GATK and dbSNP tools that produce INFO tags
that violate the VCF specification:

- `dbSNP138_ID`, `dbSNP142_ID` — present in data rows but missing from
  the `##INFO` header lines
- `1000G_phase1_release_v3_AF` — tag name starts with a digit, which is
  illegal in VCF format

These cause `bcftools norm` to fail with:
```
[E::vcf_format] Invalid BCF, the INFO tag id=38 is too large
```

The `normalize_vcf` rule handles this by running `bcftools annotate`
**first** in the pipe to strip all INFO tags before `bcftools norm` sees
them. The `--force` flag is required because the tags are missing from the
header. Using `--output-type v` (plain text VCF) at this step is also
required — BCF binary output triggers strict header validation even when
the offending tags are being removed.

The full pipe order is:
```
bcftools annotate (-O v, strip INFO) →
bcftools view (-O u, PASS + SNPs/indels + --targets) →
bcftools norm (-O u, decompose + left-align) →
bcftools filter (-O u, GQ/DP) →
bcftools view (-O z, write to disk)
```

Note: `--targets` is used instead of `--regions` for the chromosome filter
because `--regions` requires an index and cannot read from stdin in a pipe.

### Decoy Contigs (EBV / NC_007605)

PPMI samples were processed with a reference genome that includes the EBV
genome (NC_007605) and other hs37d5 decoy sequences. Using a plain hg19
reference without decoys causes `bcftools norm` to fail with:
```
[E::faidx_adjust_position] The sequence "NC_007605" was not found
```

Always use the `hs37d5` reference (see Setup Step 2). The `normalize_vcf`
rule also filters to autosomes 1–22 via `--targets` which excludes any
remaining decoy contigs from downstream analysis.

### Duplicate rsIDs After Multiallelic Split

`bcftools norm -m -both` splits multiallelic variants into biallelic
records. When both split components share the same rsID (e.g. both alleles
of a multiallelic site have `rs12564852`), PLINK2's merge step fails with:
```
Error: The biallelic variants with ID 'rs12564852' ... appear to be the
components of a 'split' multiallelic variant
```

The `merge_plink2` rule handles this with `--set-all-var-ids '@:#:$r:$a'`
which replaces all variant IDs with a unique `chr:pos:ref:alt` format
before merging.

---

## Troubleshooting

### All samples filtered by PLINK2 `--mind`

The most common failure mode. `--mind 0.10` removes samples with > 10%
missing genotypes. Diagnose first without any filter:

```bash
plink2 \
    --vcf results/merged/merged.filtered.vcf.gz \
    --missing \
    --out results/qc/debug_missing \
    --threads 4

# Distribution of per-sample missingness
awk 'NR>1 {print $5}' results/qc/debug_missing.smiss \
    | sort -n \
    | awk '{b=int($1*10)/10; c[b]++} END {for(k in c) printf "%.1f\t%d\n",k,c[k]}' \
    | sort -n
```

If most samples exceed 0.10, the issue is usually that GQ/DP filters in
`normalize_vcf` were too aggressive, leaving too many missing genotypes
per sample. Try relaxing in `config.yaml`:

```yaml
filters:
  min_gq: 10    # reduce from 20
  min_dp: 5     # reduce from 10
```

### Normalisation produces empty VCFs

Check the actual GQ/DP distribution in your data before setting thresholds:

```bash
bcftools query -f '[%GQ\n]' \
    Downloads/PPMISI001.realigned.recalibrated.g.genotyped.filtered.annotated.vcf.gz \
    | sort -n | uniq -c | head -30
```

### HWE filter removing too many variants

PPMI mixes PD cases and controls. Hardy-Weinberg equilibrium is only valid
in healthy controls — applying it to a mixed cohort removes real
disease-associated variants. Options:

```yaml
# Option A: relax threshold in config.yaml
plink:
  hwe: 1e-10    # instead of 1e-6

# Option B: apply to controls only — add to vcf_to_plink2 rule in Snakefile:
#   --hwe 1e-6 keep-fewhet
```

### Merge step runs out of memory

Increase `mem_gb_merge` in `config.yaml` and reduce `threads_merge`
(bcftools merge allocates per-thread buffers):

```yaml
resources:
  threads_merge: 4     # reduce from 8
  mem_gb_merge: 48     # increase if 32 GB is insufficient
```

### Corrupted input VCFs

Check all input VCFs for BGZF integrity before running:

```bash
for f in Downloads/PPMISI*.vcf.gz; do
    bcftools view -h "$f" > /dev/null 2>&1 || echo "CORRUPTED: $f"
done
```

Re-download any corrupted files from the PPMI data portal.

---

## Tool Versions

| Tool      | Version  | Purpose                           |
|-----------|----------|-----------------------------------|
| bcftools  | 1.19     | VCF normalisation, merge, filter  |
| htslib    | 1.19     | tabix indexing                    |
| samtools  | 1.19     | reference FASTA indexing          |
| plink2    | 2.00a5   | genotype QC and format conversion |
| snakemake | ≥ 7.0    | workflow management               |

Tool environments are defined in `envs/bcftools.yaml` and `envs/plink2.yaml`
and are created automatically by Snakemake with `--use-conda`.
