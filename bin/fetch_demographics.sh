#!/bin/bash

if [ $# -lt 1 ]
then
  echo "Usage: (peg)[1|2]"
  exit 1
fi

peg=$1
if [ $peg -eq 1 ]
then
  samplesize=580
  eth_field='Ethnicity'
elif [ $peg -eq 2 ]
then
  samplesize=209
  eth_field='EthnicityOriginal'
fi
echo "PEG ${peg}"

sql_pd_qtl << END
create temporary table tempsubjects(subject_id varchar(255),primary key(subject_id));
load data infile "/var/analysis/PD-QTLs/rawdata/merge/peg${peg}/subjectlist.txt" ignore into table tempsubjects fields terminated by '\t';

select ${peg} as study,PDstudyParkinsonsDisease as affection,Female as is_female,count(*) from peg${peg}_covariates as a,ethnicity as a2,peg_id_gwas_id_mapping as b,gwas_subject_genotypes_${samplesize} as c,datMethPEG${peg}t_sample as d where a2.peg_id=a.peg_id and b.peg_id=a.peg_id and b.dup=0 and c.gwas_id=b.gwas_id and d.sample_id=a.sample_id group by study, affection,is_female;

select ${peg} as study,PDstudyParkinsonsDisease as affection,${eth_field},count(*) from peg${peg}_covariates as a,ethnicity as a2,peg_id_gwas_id_mapping as b,gwas_subject_genotypes_${samplesize} as c,datMethPEG${peg}t_sample as d where a2.peg_id=a.peg_id and b.peg_id=a.peg_id and b.dup=0 and c.gwas_id=b.gwas_id and d.sample_id=a.sample_id group by study, affection,${eth_field};

select 'controls',avg(age),min(age),max(age) from peg${peg}_covariates as a,ethnicity as a2,peg_id_gwas_id_mapping as b,gwas_subject_genotypes_${samplesize} as c,datMethPEG${peg}t_sample as d where PDstudyParkinsonsDisease=0 and a2.peg_id=a.peg_id and b.peg_id=a.peg_id and b.dup=0 and c.gwas_id=b.gwas_id and d.sample_id=a.sample_id;
select 'cases',avg(age),min(age),max(age) from peg${peg}_covariates as a,ethnicity as a2,peg_id_gwas_id_mapping as b,gwas_subject_genotypes_${samplesize} as c,datMethPEG${peg}t_sample as d where PDstudyParkinsonsDisease=1 and a2.peg_id=a.peg_id and b.peg_id=a.peg_id and b.dup=0 and c.gwas_id=b.gwas_id and d.sample_id=a.sample_id;

select race_gwas,format(avg(k1a),3),format(avg(k2a),3),format(avg(k3a),3),format(avg(k4a),3) from peg${peg}_covariates as a,ethnicity as a2,peg_id_gwas_id_mapping as b,gwas_subject_genotypes_${samplesize} as c,datMethPEG${peg}t_sample as d where a2.peg_id=a.peg_id and b.peg_id=a.peg_id and b.dup=0 and c.gwas_id=b.gwas_id and d.sample_id=a.sample_id group by race_gwas;
select race_q,format(avg(k1a),3),format(avg(k2a),3),format(avg(k3a),3),format(avg(k4a),3) from peg${peg}_covariates as a,ethnicity as a2,peg_id_gwas_id_mapping as b,gwas_subject_genotypes_${samplesize} as c,datMethPEG${peg}t_sample as d where a2.peg_id=a.peg_id and b.peg_id=a.peg_id and b.dup=0 and c.gwas_id=b.gwas_id and d.sample_id=a.sample_id group by race_q;

#select race_gwas,format(avg(k1a),3),format(avg(k2a),3),format(avg(k3a),3),format(avg(k4a),3) from tempsubjects as a,ethnicity as b where b.peg_id=a.subject_id group by race_gwas;
#select race_q,format(avg(k1a),3),format(avg(k2a),3),format(avg(k3a),3),format(avg(k4a),3) from tempsubjects as a,ethnicity as b where b.peg_id=a.subject_id group by race_q;
END
