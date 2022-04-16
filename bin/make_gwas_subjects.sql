use pd_qtl;
drop table if exists  gwas_subjects ;
create table  gwas_subjects  (seq int auto_increment, gwas_id  varchar(255),  primary key ( seq ),unique key index_gwas_id(gwas_id));
load data infile '/home/garyc/analysis/PD-QTLs/rawdata/Genetics/subjects.txt' into table  gwas_subjects  fields terminated by '\t' ignore 0 lines(gwas_id);
