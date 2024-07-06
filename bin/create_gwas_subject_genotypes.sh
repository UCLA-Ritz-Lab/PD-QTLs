#!/bin/bash

if [ $# -lt 1 ]
then
  echo "Usage: [peg1cases|peg1controls|peg2cases]"
  exit 1
fi
samplesize=$1

sql_pd_qtl << END
use pd_qtl;
drop table if exists gwas_subject_genotypes_${samplesize} ;
create table  gwas_subject_genotypes_${samplesize}  (gwas_id varchar(25), genotype_string mediumtext, primary key(gwas_id));
load data infile "/home/garyc/analysis/PD-QTLs/rawdata/Genetics/gwas_additive_genotypes_${samplesize}.txt" into table  gwas_subject_genotypes_${samplesize}  fields terminated by '\t' ignore 1 lines;
END
exit 0
