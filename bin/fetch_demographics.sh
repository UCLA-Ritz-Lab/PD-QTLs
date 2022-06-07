#!/bin/bash

if [ $# -lt 1 ]
then
  echo "Usage: (peg)[1|2]"
  exit 1
fi

peg=$1
sql_pd_qtl << END
select "PEG ${peg}";
create temporary table tempsubjects(subject_id varchar(255),primary key(subject_id));
load data infile "/var/analysis/PD-QTLs/rawdata/merge/peg${peg}/subjectlist.txt" ignore into table tempsubjects fields terminated by '\t';
select "age";
select min(age),avg(age),max(age) from tempsubjects as a,peg${peg}_covariates as b where b.peg_id=a.subject_id;
select "female";
select female,count(*) from tempsubjects as a,peg${peg}_covariates as b where b.peg_id=a.subject_id group by female;
END
