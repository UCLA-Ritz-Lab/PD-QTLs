use pd_qtl;
drop table if exists  genotypes ;
create table  genotypes  (snpid varchar(255), genotype_string mediumtext, primary key(snpid));
load data infile '/home/garyc/analysis/PD-QTLs/rawdata/Genetics/genotypes.txt' into table  genotypes  fields terminated by '\t' ignore 0 lines;
