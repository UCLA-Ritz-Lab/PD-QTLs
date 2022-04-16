use pd_qtl;
drop table if exists  sampleid_gwasid_mapping ;
create table  sampleid_gwasid_mapping  (peg_id  varchar(25), gwas_id  varchar(25), dup  tinyint,  key ( peg_id ),  key index_gwas_id(gwas_id));
load data infile '/home/garyc/analysis/PD-QTLs/rawdata/Genetics/sampleid_gwasid_mapping_import.csv' into table sampleid_gwasid_mapping  fields terminated by ',' ignore 1 lines;
