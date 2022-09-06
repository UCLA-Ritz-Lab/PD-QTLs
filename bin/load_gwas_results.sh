#!/bin/bash

if [ $# -lt 1 ]
then
  echo "Usage: <plink2 output>"
  exit 1
fi
plinkout=$1

sql_pd_qtl << END
drop table if exists peg1_gwas;
create table peg1_gwas(CHROM varchar(30),POS int unsigned,snp_id varchar(255), ref varchar(10),  alt varchar(10), allele1 varchar(10), firth_status varchar(1), test varchar(255),  obs_ct smallint unsigned,  odd_ratio float,beta_se float, z_stat float, p_value float, error_code varchar(10),primary key(snp_id,test));
load data infile "/var/analysis/PD-QTLs/rawdata/Genetics/${plinkout}" into table peg1_gwas fields terminated by '\t' ignore 1 lines;
END
