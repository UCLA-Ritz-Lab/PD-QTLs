use pd_qtl;
drop table if exists  peg2_covariates ;
create table  peg2_covariates  (sample_id  varchar(25), peg_id varchar(25)  , Female tinyint unsigned, Age  float, PlasmaBlast  float, CD8pCD28nCD45RAn  float, CD8_Naive  float, CD4_Naive  float, CD8T  float, CD4T  float, NK  float, Bcell  float, Mono  float, Gran  float, EthnicityAIM_RF text, RFvoteAsian text, RFvoteBlack  text, RFvoteCaucasian  text, RFvoteHispanic  text, EthnicityOriginal  text,  primary key ( sample_id ), unique key (peg_id));
load data infile '/home/garyc/analysis/PD-QTLs/rawdata/Methylation/peg2_methylation_controlvar_import.csv' into table  peg2_covariates  fields terminated by ',' ignore 1 lines;
alter table peg2_covariates add column PDstudyParkinsonsDisease tinyint unsigned;
create temporary table disease_status(peg_id varchar(25),disease_status tinyint unsigned,primary key(peg_id));
load data infile '/home/garyc/analysis/PD-QTLs/rawdata/peg_disease_status.csv'  into table disease_status fields terminated by ',' ignore 1 lines;
update peg2_covariates as a,disease_status as b set a.PDstudyParkinsonsDisease=b.disease_status where a.peg_id=b.peg_id;
