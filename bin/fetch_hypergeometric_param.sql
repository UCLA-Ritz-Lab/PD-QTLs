create temporary table all_snps(chrom varchar(10),snpid varchar(30),a tinyint unsigned, position int unsigned, allele1 varchar(10),allele2 varchar(10),primary key(snpid));
load data infile '/var/analysis/PD-QTLs/rawdata/Genetics/PEG.phased.580.methex.bim' into table all_snps fields terminated by '\t';

create temporary table me_qtls(eqtl_type varchar(30),snpid varchar(30),effect_allele varchar(30), probe varchar(30),statistic float, pvalue float, fdr float, beta float,primary key(eqtl_type,snpid,probe));
load data infile '/var/analysis/PD-QTLs/rawdata/merge/cis_eqtls.txt' into table me_qtls fields terminated by '\t';
load data infile '/var/analysis/PD-QTLs/rawdata/merge/trans_eqtls.txt' into table me_qtls fields terminated by '\t';
load data infile '/var/analysis/PD-QTLs/rawdata/merge/all_eqtls.txt' into table me_qtls fields terminated by '\t';

create temporary table  pd_snps_nall  (snpid  varchar(255), se float, beta float, effect_allele varchar(10), other_allele varchar(10), eaf float, phenotype text,chrom varchar(30), position int unsigned, samplesize mediumint unsigned, ncase mediumint unsigned, ncontrols mediumint unsigned, pvalue float, units text, gene text,  primary key ( snpid ));
load data infile '/home/garyc/analysis/PD-QTLs/rawdata/Genetics/nalls_pd_snps.txt' into table pd_snps_nall fields terminated by '\t' ignore 1 lines;

select "ALL SNPS";
select count(snpid) from all_snps;
select "ALL cis MEQTLS";
select count(distinct b.snpid) from all_snps as a,me_qtls as b where b.snpid=a.snpid and b.eqtl_type='cis';
select "ALL PD SNPS";
select count(a.snpid) from pd_snps_nall as a,all_snps as b where b.snpid=a.snpid;
select "ALL PD SNPS that are cis MEQTLS";
select count(distinct b.snpid) from all_snps as a,me_qtls as b,pd_snps_nall as c where b.snpid=a.snpid and b.eqtl_type='cis' and c.snpid=a.snpid;

select "ALL SNPS";
select count(snpid) from all_snps;
select "ALL trans MEQTLS";
select count(distinct b.snpid) from all_snps as a,me_qtls as b where b.snpid=a.snpid and b.eqtl_type='trans';
select "ALL PD SNPS";
select count(a.snpid) from pd_snps_nall as a,all_snps as b where b.snpid=a.snpid;
select "ALL PD SNPS that are trans MEQTLS";
select count(distinct b.snpid) from all_snps as a,me_qtls as b,pd_snps_nall as c where b.snpid=a.snpid and b.eqtl_type='trans' and c.snpid=a.snpid;

select "ALL SNPS";
select count(snpid) from all_snps;
select "ALL cistrans MEQTLS";
select count(distinct b.snpid) from all_snps as a,me_qtls as b where b.snpid=a.snpid and b.eqtl_type='all';
select "ALL PD SNPS";
select count(a.snpid) from pd_snps_nall as a,all_snps as b where b.snpid=a.snpid;
select "ALL PD SNPS that are cistrans MEQTLS";
select count(distinct b.snpid) from all_snps as a,me_qtls as b,pd_snps_nall as c where b.snpid=a.snpid and b.eqtl_type='all' and c.snpid=a.snpid;
