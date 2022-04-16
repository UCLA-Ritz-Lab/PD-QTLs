use pd_qtl;
drop table if exists  peg_id_gwas_id_mapping ;
create table  peg_id_gwas_id_mapping  (peg_id varchar(25), gwas_id varchar(25), dup tinyint,  key index_peg_id(peg_id), key index_gwas_id(gwas_id));
load data infile '/home/garyc/analysis/PD-QTLs/rawdata/Genetics/peg_id_gwas_id_mapping_import.csv' into table peg_id_gwas_id_mapping  fields terminated by ',' ignore 1 lines;
