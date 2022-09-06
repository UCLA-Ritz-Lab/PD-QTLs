drop table if exists pd_probes;
create table pd_probes(IlmnID varchar(255),primary key(IlmnID));
load data infile '/var/analysis/PD-QTLs/rawdata/Methylation/pd_probes.txt' ignore into table pd_probes fields terminated by '\t' ignore 1 lines;
