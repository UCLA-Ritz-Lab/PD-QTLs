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
select a.peg_id,concat(c18_mapping.sample_id,'_',hilic_mapping.sample_id) as c18_hilic_sample_id,Female,Age,k1a,k2a,k3a,k4a,PDstudyParkinsonsDisease,b.gwas_id,c.genotype_string,concat(c18.betas_string,',',hilic.betas_string) as c18_hilic_rawdata from peg${peg}_covariates as a,ethnicity as a2,peg_id_gwas_id_mapping as b,gwas_subject_genotypes_${samplesize} as c,peg_metabolome_rawdata as c18,peg_metabolome_rawdata as hilic,peg_metabolome_sample_mapping as c18_mapping, peg_metabolome_sample_mapping as hilic_mapping where a.PDstudyParkinsonsDisease=${disease_status} and a2.peg_id=a.peg_id and b.peg_id=a.peg_id and b.dup=0 and c.gwas_id=b.gwas_id and c18_mapping.sample_id=c18.sample_id and c18.dataset='c18' and c18_mapping.peg_id=a.peg_id and hilic_mapping.sample_id=hilic.sample_id and hilic.dataset='hilic' and hilic_mapping.peg_id=a.peg_id
END
