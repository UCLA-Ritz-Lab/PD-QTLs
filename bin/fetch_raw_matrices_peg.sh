#!/bin/bash

if [ $# -lt 3 ]
then
  echo "Usage: [1|2 (PEG1 or PEG2)] [peg1cases|peg1controls|peg2cases] [0|1 (disease status)]"
  exit 1
fi
peg=$1
samplesize=$2
disease_status=$3

sql_pd_qtl << END
select a.peg_id,a.sample_id,Female,Age,k1a,k2a,k3a,k4a,Mono,Gran,CD4T,NK,CD8_naive,CD8pCD28nCD45RAn,PlasmaBlast,b.gwas_id,c.genotype_string,d.betas_string from peg${peg}_covariates as a,ethnicity as a2,peg_id_gwas_id_mapping as b,gwas_subject_genotypes_${samplesize} as c,datMethPEG${peg}t_sample as d where a.PDstudyParkinsonsDisease=${disease_status} and a2.peg_id=a.peg_id and b.peg_id=a.peg_id and b.dup=0 and c.gwas_id=b.gwas_id and d.sample_id=a.sample_id;
END
