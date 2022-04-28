use pd_qtl;
drop table if exists  illumina_annotation ;
create table  illumina_annotation  (IlmnID  varchar(255), chrom varchar(25),mapinfo int unsigned, probe_snps varchar(50), probe_snps_10 varchar(50),  primary key ( IlmnID ),unique key index_position(chrom,mapinfo));
load data infile '/home/garyc/analysis/PD-QTLs/rawdata/Methylation/illumina_annotation_import.csv' into table  illumina_annotation  fields terminated by ',' ignore 1 lines;
