#!/bin/bash

tsv_files=`ls dat*tsv`
for tsv_file in ${tsv_files}
do
  keytype=`echo $tsv_file|cut -f2 -d'_'|cut -f1 -d\.`
  echo Running "../../bin/tsv2sql.sh ${tsv_file} ${keytype}_id"
  ../../bin/tsv2sql.sh ${tsv_file} ${keytype}_id
done
