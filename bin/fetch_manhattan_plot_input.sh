#!/bin/bash

if [ $# -lt  2 ]
then
  echo "Usage: <peg study (peg1|peg2)> <me_qtl_type (cis|trans)>" 
  exit 1
fi

peg=$1
me_qtl_type=$2

sql_pd_qtl << END
select chrom as CHR,pos as BP,snpid as SNP,pvalue as P from snpinfo as a,me_qtls as b where b.snp_id=a.snpid and peg="${peg}" and me_qtl_type="${me_qtl_type}" order by pvalue;
END

