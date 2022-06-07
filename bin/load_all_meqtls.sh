#!/bin/bash

sql_pd_qtl  << END
drop table if exists me_qtls;
create table me_qtls(peg varchar(10),me_qtl_type varchar(10),snp_id varchar(55), allele varchar(255),gene varchar(55),statistic float,pvalue double,FDR double,beta float,primary key(peg,me_qtl_type,snp_id,allele,gene), index pairing(snp_id,allele,gene));
END
pegs='peg1 peg2'
meqtl_types='cis trans'
for peg in $pegs
do
  for meqtl_type in $meqtl_types
  do
    echo "$peg $meqtl_type"
    sed "s/^/${peg}\t/" /var/analysis/PD-QTLs/rawdata/merge/${peg}/${meqtl_type}_eqtls.txt > /var/analysis/PD-QTLs/rawdata/import.tmp
    echo "load data infile '/var/analysis/PD-QTLs/rawdata/import.tmp' into table me_qtls;" | sql_pd_qtl
  done
done
