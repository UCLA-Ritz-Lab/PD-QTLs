#!/bin/bash

project_dir='/home/garyc/analysis/PD-QTLs'
bindir=${project_dir}/bin

if [ $# -lt 2 ]
then
  echo "<tsv_filename> <unique_col>"
  exit 1
fi

tsv_filename=$1
unique_col=$2
tablename=`echo $tsv_filename|cut -d\. -f1`

echo Importing into table $tablename

string='databases';

#echo "create table ${tablename}(seq int auto_increment,${unique_col} varchar(50), csv_string text, primary key(seq),unique key index_${unique_col}(${unique_col}));"
sql_pd_qtl << END
drop table if exists ${tablename};
#create table ${tablename}(seq int auto_increment,subject_id varchar(50), csv_string text, primary key(seq));
create table ${tablename}(seq int auto_increment,${unique_col} varchar(50), csv_string mediumtext, primary key(seq),unique key index_${unique_col}(${unique_col}));
load data infile "${PWD}/${tsv_filename}" into table ${tablename} fields terminated by '\t'(${unique_col},csv_string);
END
#echo Converting ${robject}
#sed "s/ROBJECT/${robject}/g" ${bindir}/import_rdata.template | sed "s/NOTRANSPOSE/${notranspose}/g" | sed "s/WITHTRANSPOSE/${withtranspose}/g" > import_rdata.r
#R --no-save< import_rdata.r
#sed "s/,/\t/" tmp_${notranspose}.csv > ${robject}_${notranspose}.tsv
#sed "s/,/\t/" tmp_${withtranspose}.csv > ${robject}_${withtranspose}.tsv
exit 0
