#!/bin/bash

if [ $# -lt 1 ]
then
  echo "Usage: [cis|trans|all]"
  exit 1
fi

eqtl_type=$1
sql_pd_qtl << END
#create temporary table pd_probes(IlmnID varchar(255),primary key(IlmnID));
#load data infile '/var/analysis/PD-QTLs/rawdata/Methylation/pd_probes.txt' ignore into table pd_probes fields terminated by '\t' ignore 1 lines;

#create temporary table all_snps(chrom varchar(10),snpid varchar(30),a tinyint unsigned, position int unsigned, allele1 varchar(10),allele2 varchar(10),primary key(snpid));
#load data infile '/var/analysis/PD-QTLs/rawdata/Genetics/PEG.phased.580.methex.bim' into table all_snps fields terminated by '\t';

create temporary table me_qtls_peg1(eqtl_type varchar(30),snpid varchar(255),effect_allele varchar(30), probe varchar(30),statistic float, pvalue float, fdr float, beta float,primary key(eqtl_type,snpid,effect_allele,probe));
load data infile "/var/analysis/PD-QTLs/rawdata/merge/peg1/${eqtl_type}_eqtls.txt" into table me_qtls_peg1 fields terminated by '\t';

create temporary table me_qtls_peg2(eqtl_type varchar(30),snpid varchar(255),effect_allele varchar(30), probe varchar(30),statistic float, pvalue float, fdr float, beta float,primary key(eqtl_type,snpid,effect_allele,probe));
load data infile "/var/analysis/PD-QTLs/rawdata/merge/peg2/${eqtl_type}_eqtls.txt" into table me_qtls_peg2 fields terminated by '\t';

select * from me_qtls_peg1 as a,me_qtls_peg2 as b where a.eqtl_type=b.eqtl_type and a.snpid=b.snpid and a.effect_allele=b.effect_allele and a.probe=b.probe;
END



