create temporary table kimberly(gene varchar(55),primary key(gene));
load data infile '/var/analysis/PD-QTLs/rawdata/Methylation/kimberly_30_cpg.txt' into table kimberly;
select * from me_qtls as a,kimberly as b where a.gene=b.gene;

