#!/bin/bash

if [ $# -lt 1 ]
then
  echo "Usage: [1|2] (peg1|2)"
  exit 1
fi

peg=$1

sql_pd_qtl << END
drop table if exists cis_eqtls;
create table cis_eqtls(eqtl_type varchar(30),snpid varchar(255),ref_allele varchar(50),gene varchar(255),statistic float, pvalue float, fdr float, beta float, primary key(eqtl_type,snpid,ref_allele,gene));
load data infile "/home/garyc/analysis/PD-QTLs/rawdata/merge/peg${peg}/cis_eqtls.txt" into table cis_eqtls fields terminated by '\t';
select chrom,pos,pos from snpinfo as a,cis_eqtls as b where b.eqtl_type='cis' and b.snpid=a.snpid group by chrom,pos order by chrom,pos;
END
