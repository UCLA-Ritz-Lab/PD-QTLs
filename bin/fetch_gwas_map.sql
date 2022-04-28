#create temporary table gwas_snp(snpid varchar(255),primary key(snpid));
#load data infile '/home/garyc/analysis/PD-QTLs/rawdata/merge/snplist_t.txt' into table gwas_snp fields terminated by ',';
select a.snpid as snp, ref_allele,chrom as chr, pos from snpinfo as a;
#select a.snpid as snp, ref_allele,chrom as chr, pos from snpinfo as a,gwas_snp as b where b.snpid=a.snpid;
#drop table gwas_snp;

