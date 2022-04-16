use pd_qtl;
drop table if exists  gwas_subject_genotypes ;
create table  gwas_subject_genotypes  (gwas_id varchar(25), genotype_string mediumtext, primary key(gwas_id));
load data infile '/home/garyc/analysis/PD-QTLs/rawdata/Genetics/gwas_additive_genotypes.txt' into table  gwas_subject_genotypes  fields terminated by '\t' ignore 1 lines;
