#!/bin/bash

if [ $# -lt 1 ]
then
  echo "Usage: [tablename: me_qtls|m_qtls] [qtl_type: all|cis|trans] [peg: peg1|peg2|all]"
  exit 1
fi

tablename=$1
qtl_type=$2
peg=$3

if [ $peg == 'all' ]
then
  sql_string="select chrom,pos,pos from snpinfo as a,${tablename} as b1,${tablename} as b2 where b1.peg='peg1' and b1.qtl_type='${qtl_type}' and b1.snp_id=a.snpid and b2.peg='peg2' and b2.qtl_type=b1.qtl_type and b2.snp_id=b1.snp_id and b2.allele=b1.allele and b2.gene=b1.gene group by chrom,pos order by chrom,pos;"
else
  sql_string="select chrom,pos,pos from snpinfo as a,${tablename} as b where b.peg='${peg}' and b.qtl_type='${qtl_type}' and b.snp_id=a.snpid group by chrom,pos order by chrom,pos;"
fi

echo $sql_string
echo $sql_string | sql_pd_qtl
