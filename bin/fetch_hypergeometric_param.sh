#!/bin/bash

if [ $# -lt 3 ]
then
  echo "Usage: [qtl tablename: me_qtls|m_qtls] [qtl_type: cis|trans|all]  [peg: peg1cases|peg1controls|peg2cases]"
  exit 1
fi

qtl_tablename=$1
qtl_type=$2
peg=$3
sql_pd_qtl << END
create temporary table all_snps_${peg}(snpid varchar(255),primary key(snpid));
load data infile "/var/analysis/PD-QTLs/rawdata/merge/${peg}/snplist_t.txt" into table all_snps_${peg} fields terminated by '\t';

#select "GWAS SNPs";
select count(snpid) as gwas_snps from all_snps_${peg};
#select "GWAS SNPs that are ${qtl_type} QTLS";
select count(distinct b.snp_id) as gwas_eqtls from all_snps_${peg} as a,${qtl_tablename} as b where b.peg="${peg}" and b.snp_id=a.snpid and b.qtl_type="${qtl_type}";
#select "PD risk SNPS";
select count(a.snpid) as pd_snps from pd_snps_nall as a,all_snps_${peg} as b,snpinfo as c where a.snpid=c.snpid and c.snpid2=b.snpid;
#select "PD risk SNPS that are ${qtl_type} QTLS";
select count(distinct b.snp_id) as pd_qtls from all_snps_${peg} as a,${qtl_tablename} as b,pd_snps_nall as c,snpinfo as d where b.peg="${peg}" and  b.snp_id=a.snpid and b.qtl_type="${qtl_type}" and c.snpid=d.snpid and  d.snpid2=a.snpid;
END
