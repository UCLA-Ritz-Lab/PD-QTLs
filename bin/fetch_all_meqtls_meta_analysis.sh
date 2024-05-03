#!/bin/bash

if [ $# -lt 2 ]
then
  echo "Usage: [peg1|peg2] [cis|trans]"
  exit 1
fi
peg=$1
cistrans=$2

sql_pd_qtl << END
select concat(snp_id,',',allele,',',gene) as 'snp,allele,gene',ref_allele,alt_allele,beta,beta/statistic as se,pvalue from me_qtls as a, snpinfo as b where b.snpid=a.snp_id and a.peg='${peg}' and a.qtl_type='${cistrans}' 
END
