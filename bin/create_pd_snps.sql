drop table if exists pd_snps_nall;
create table  pd_snps_nall  (snpid  varchar(255), se float, beta float, effect_allele varchar(10), other_allele varchar(10), eaf float, phenotype text,chrom varchar(30), position int unsigned, samplesize mediumint unsigned, ncase mediumint unsigned, ncontrols mediumint unsigned, pvalue float, units text, gene text,  primary key ( snpid ));
load data infile '/home/garyc/analysis/PD-QTLs/rawdata/Genetics/nalls_pd_snps.txt' into table pd_snps_nall fields terminated by '\t' ignore 1 lines;
