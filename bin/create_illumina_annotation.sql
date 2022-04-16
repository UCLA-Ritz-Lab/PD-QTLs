use pd_qtl;
drop table if exists  illumina_annotation ;
create table  illumina_annotation  (IlmnID  varchar(255), chrom varchar(25), mapinfo int unsigned,  primary key ( IlmnID ), unique key index_position(chrom,mapinfo));
load data infile '/home/garyc/analysis/PD-QTLs/rawdata/Methylation/illumina_annotation_import.csv' into table  illumina_annotation  fields terminated by ',' ignore 1 lines;
