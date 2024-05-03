#!/bin/bash

if [ $# -lt  2 ]
then
  echo "Usage: <peg study (peg1|peg2|all)> <me_qtl_type (cis|trans)>" 
  exit 1
fi

peg=$1
me_qtl_type=$2

if [ $peg == 'all' ]
then
  sql_string="select chrom as CHR,pos as BP,snpid as SNP,c.pvalue as P,b1.gene,b2.allele,b2.* from snpinfo as a,me_qtls as b1,me_qtls as b2,meta_analysis as c where b1.snp_id=a.snpid and b1.peg='peg1' and b1.qtl_type='${me_qtl_type}' and b1.snp_id=a.snpid and b2.peg='peg2' and b2.qtl_type=b1.qtl_type and b2.snp_id=b1.snp_id and b1.allele=b2.allele and b1.gene=b2.gene and c.gene=b1.gene and c.qtl_type=b1.qtl_type and c.snp_id=a.snpid order by c.pvalue limit 10;"
else
  sql_string="select chrom as CHR,pos as BP,snpid as SNP,c.pvalue as P from snpinfo as a,me_qtls as b1,meta_analysis as c where b1.snp_id=a.snpid and b1.peg='${peg}' and b1.qtl_type='${me_qtl_type}' and c.qtl_type='${me_qtl_type}' and c.snp_id=a.snpid order by c.pvalue;"
fi


echo $sql_string  
echo $sql_string |sql_pd_qtl
