use pd_qtl;
drop table if exists ethnicity;
create table ethnicity (peg_id varchar(25),race_gwas varchar(255),race_q varchar(255), k1a float, k2a float, k3a float, k4a float ,primary key(peg_id));
load data infile '/home/garyc/analysis/PD-QTLs/rawdata/Ethnicity/gc_import_revised.csv' into table ethnicity fields terminated by ',' ignore 1 lines;
