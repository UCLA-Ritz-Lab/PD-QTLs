#!/bin/bash

sql_pd_qtl << END
drop table if exists qtl_hotspots;
create table qtl_hotspots(qtl_type varchar(20),peg varchar(10),qtl varchar(255),total int unsigned,primary key(qtl_type,peg,qtl),key index_qtl(qtl));
load data infile '/var/analysis/PD-QTLs/rawdata/merge/peg1/trans_hotspots.txt' into table qtl_hotspots;
load data infile '/var/analysis/PD-QTLs/rawdata/merge/peg2/trans_hotspots.txt' into table qtl_hotspots;
load data infile '/var/analysis/PD-QTLs/rawdata/merge_metabolome/peg1/metab_hotspots.txt' into table qtl_hotspots;
load data infile '/var/analysis/PD-QTLs/rawdata/merge_metabolome/peg2/metab_hotspots.txt' into table qtl_hotspots;
END
