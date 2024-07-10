#!/bin/bash

if [ $# -lt  2 ]
then
  echo "Usage: <peg study (peg1|peg2|all)> <me_qtl_type (cis|trans)>" 
  exit 1
fi

peg=$1
me_qtl_type=$2

sql_string="select chrom as CHR,pos as BP,snpid as SNP,c.pvalue as P from snpinfo as a,me_qtls as b1,meta_analysis as c where b1.snp_id=a.snpid2 and b1.peg='${peg}' and b1.qtl_type='${me_qtl_type}' and c.qtl_type='${me_qtl_type}' and c.snp_id=a.snpid2 order by c.pvalue;"


#echo $sql_string  
echo $sql_string |sql_pd_qtl
