#!/bin/bash

if [ $# -lt 1 ]
then
  echo "Usage: (peg)[1|2]"
  exit 1
fi

peg=$1
echo "PEG ${peg}"
sql_pd_qtl << END
create temporary table tempsubjects(subject_id varchar(255),primary key(subject_id));
load data infile "/var/analysis/PD-QTLs/rawdata/merge/peg${peg}/subjectlist.txt" ignore into table tempsubjects fields terminated by '\t';
select format(min(age),3),format(avg(age),3),format(max(age),3) from tempsubjects as a,peg${peg}_covariates as b where b.peg_id=a.subject_id;
select female,count(*) from tempsubjects as a,peg${peg}_covariates as b where b.peg_id=a.subject_id group by female;
select race_gwas,format(avg(k1a),3),format(avg(k2a),3),format(avg(k3a),3),format(avg(k4a),3) from tempsubjects as a,ethnicity as b where b.peg_id=a.subject_id group by race_gwas;
select race_q,format(avg(k1a),3),format(avg(k2a),3),format(avg(k3a),3),format(avg(k4a),3) from tempsubjects as a,ethnicity as b where b.peg_id=a.subject_id group by race_q;
END
