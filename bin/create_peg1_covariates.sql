use pd_qtl;
drop table if exists  peg1_covariates ;
create table  peg1_covariates  (peg_id  varchar(25), sample_id  varchar(25), Female tinyint unsigned, Age  float, Ethnicity  text, PlasmaBlast  float, CD8pCD28nCD45RAn  float, CD8_naive  float, CD4T  float, NK  float, Mono  float, Gran  float, PDstudyParkinsonsDisease  tinyint, RFvoteHispanic  float,  primary key (peg_id), unique key index_sample_id(sample_id));
load data infile '/home/garyc/analysis/PD-QTLs/rawdata/Methylation/peg1_import.csv' into table  peg1_covariates  fields terminated by ',' ignore 1 lines;
