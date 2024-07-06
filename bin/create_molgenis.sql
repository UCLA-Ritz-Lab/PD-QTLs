use pd_qtl;
drop table if exists  molgenis;
create table  molgenis  (pvalue float,snpid varchar(255),probename varchar(255),cistrans varchar(10),effect_allele varchar(10),zscore float,fdr float,primary key(snpid,probename));
load data infile '/home/garyc/analysis/PD-QTLs/rawdata/external_qtl/molgenis_import.txt' into table molgenis fields terminated by '\t' ignore 1 lines;
