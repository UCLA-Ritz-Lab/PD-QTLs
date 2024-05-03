#!/bin/bash

sql_pd_qtl  << END
drop table if exists meta_analysis;
create table meta_analysis(qtl_type varchar(10),snp_id varchar(55), gene varchar(55),zscore float,pvalue double,primary key(qtl_type,snp_id,gene));
END
echo "load data infile '/var/analysis/PD-QTLs/rawdata/merge/metal_import.txt' into table meta_analysis;" | sql_pd_qtl
