#!/bin/bash

sql_pd_qtl  << END
drop table if exists met_qtls;
create table met_qtls(peg varchar(20),qtl_type varchar(10),snp_id varchar(55), gene varchar(55),statistic float,pvalue double,FDR double,beta float,primary key(peg,qtl_type,snp_id,gene), index pairing(snp_id,gene));
END
pegs='peg1cases peg1controls peg2cases'
mqtl_types='all'
for peg in $pegs
do
  for mqtl_type in $mqtl_types
  do
    echo "$peg $mqtl_type"
    sed "s/^/${peg}\t/" /var/analysis/PD-QTLs/rawdata/merge_metabolome/${peg}/${mqtl_type}_eqtls.txt > /var/analysis/PD-QTLs/rawdata/import.tmp
    echo "load data infile '/var/analysis/PD-QTLs/rawdata/import.tmp' into table met_qtls;" | sql_pd_qtl
  done
done
