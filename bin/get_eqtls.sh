#!/bin/bash

if [ $# -lt 1 ]
then
  echo "Usage: [tablename: me_qtls|m_qtls] [qtl_type: all|cis|trans] [peg: peg1|peg2]"
  exit 1
fi

tablename=$1
qtl_type=$2
peg=$3

sql_pd_qtl << END
select chrom,pos,pos from snpinfo as a,${tablename} as b where b.peg="${peg}" and b.qtl_type="${qtl_type}" and b.snp_id=a.snpid group by chrom,pos order by chrom,pos;
END
