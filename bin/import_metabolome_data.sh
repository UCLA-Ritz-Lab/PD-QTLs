#!/bin/bash

project_dir='/home/garyc/analysis/PD-QTLs'
bindir=${project_dir}/bin
mapping_table='peg_metabolome_sample_mapping'

sql_pd_qtl << END1
drop table if exists ${mapping_table};
create table ${mapping_table}(sample_id varchar(50),peg_id varchar(50),primary key(sample_id),key index_peg_id(peg_id));
load data infile "${PWD}/sample_map_c18_import.tsv" into table ${mapping_table} fields terminated by '\t'(sample_id,peg_id);
load data infile "${PWD}/sample_map_hilic_import.tsv" into table ${mapping_table} fields terminated by '\t'(sample_id,peg_id);
END1

rawdata_table='peg_metabolome_rawdata'
sql_pd_qtl << END2
drop table if exists ${rawdata_table};
create table ${rawdata_table}(seq int auto_increment,dataset varchar(20),sample_id varchar(50),betas_string mediumtext,primary key(seq),key index_dataset(dataset),unique key index_sample_id(sample_id));
load data infile "${PWD}/c18_import.tsv" into table ${rawdata_table} fields terminated by '\t'(dataset,sample_id,betas_string);
load data infile "${PWD}/hilic_import.tsv" into table ${rawdata_table} fields terminated by '\t'(dataset,sample_id,betas_string);
END2

