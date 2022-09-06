#!/bin/bash

if [ $# -lt 1 ]
then
  echo "Usage: <famfile>"
  exit 1
fi
famfile=$1

sql_pd_qtl << END
create temporary table famfile(fid varchar(20),iid varchar(20),a tinyint, b tinyint, c tinyint, affection tinyint, primary key(fid,iid));
load data infile "/var/analysis/PD-QTLs/rawdata/Genetics/${famfile}" into table famfile fields terminated by '\t';
select famfile.fid as FID,famfile.iid as IID,cov.female,cov.age,cov.RFvoteHispanic from famfile left join peg_id_gwas_id_mapping as mapping on mapping.gwas_id=famfile.iid inner join peg1_covariates as cov on cov.peg_id =  mapping.peg_id;
END
