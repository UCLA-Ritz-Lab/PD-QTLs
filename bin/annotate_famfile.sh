#!/bin/bash

if [ $# -lt 1 ]
then
  echo "Usage: <famfile>"
  exit 1
fi
famfile=$1

sql_pd_qtl << END
  create temporary table famfile(FID varchar(20),IID varchar(20),a tinyint, b tinyint, c tinyint, affection tinyint, primary key(fid,iid));
  load data infile "/var/analysis/PD-QTLs/rawdata/Genetics/${famfile}" into table famfile fields terminated by '\t';
  update famfile as fam,peg_id_gwas_id_mapping as mapping,peg1_covariates as peg1 set fam.affection=(peg1.PDstudyParkinsonsDisease+1) where mapping.peg_id=peg1.peg_id and fam.iid=mapping.gwas_id;
  select famfile.FID as FID, famfile.IID as IID,famfile.a,famfile.b,famfile.c,famfile.affection from famfile;
END
