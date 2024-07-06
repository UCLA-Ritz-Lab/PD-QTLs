use pd_qtl;
drop table if exists  snpinfo ;
create table  snpinfo  (snpid varchar(255),chrom varchar(10),pos int unsigned, ref_allele varchar(10),alt_allele varchar(10),primary key(snpid),unique key index_position(chrom,pos));
load data infile '/home/garyc/analysis/PD-QTLs/rawdata/Genetics/snpinfo.txt' into table  snpinfo  fields terminated by '\t' ignore 0 lines;
alter table snpinfo add snpid2 varchar(255);
create unique index index_snpid2 on snpinfo(snpid2);
update snpinfo set snpid2=concat(snpid,'_',ref_allele);
