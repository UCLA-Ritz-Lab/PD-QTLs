#!/bin/bash
set -euo pipefail

META_DIR=../../results/META_MQTL/meta
GENO_DIR=../../results/PPMI_MQTL/genotypes
BETA=../../results/PPMI_GCP/local_results/final/beta_matrix.csv.gz
OUT_DIR=../../results/PPMI_LONG

for GRP in cases controls; do
  echo "=== Preparing $GRP ==="
  HITS=${META_DIR}/${GRP}_sig_fdr05.csv

  tail -n +2 "$HITS" | awk -F',' '{print $2}' | sort -u > ${OUT_DIR}/${GRP}_target_snps.txt
  tail -n +2 "$HITS" | awk -F',' '{print $3}' | sort -u > ${OUT_DIR}/${GRP}_target_cpgs.txt

  echo "  target SNPs: $(wc -l < ${OUT_DIR}/${GRP}_target_snps.txt)"
  echo "  target CpGs: $(wc -l < ${OUT_DIR}/${GRP}_target_cpgs.txt)"

  zcat ${GENO_DIR}/ppmi_${GRP}.geno.tsv.gz | \
    awk -F'\t' 'NR==FNR{t[$1]=1; next} FNR==1 || ($1 in t)' ${OUT_DIR}/${GRP}_target_snps.txt - \
    > ${OUT_DIR}/${GRP}_geno_subset_full.tsv

  zcat "$BETA" | \
    awk -F',' 'NR==FNR{t[$1]=1; next} FNR==1 || ($1 in t)' ${OUT_DIR}/${GRP}_target_cpgs.txt - \
    > ${OUT_DIR}/${GRP}_beta_subset_full.csv

  echo "  geno subset rows: $(wc -l < ${OUT_DIR}/${GRP}_geno_subset_full.tsv)"
  echo "  beta subset rows: $(wc -l < ${OUT_DIR}/${GRP}_beta_subset_full.csv)"
done
