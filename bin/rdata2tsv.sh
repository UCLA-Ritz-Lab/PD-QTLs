#!/bin/bash

project_dir='/home/garyc/analysis/PD-QTLs'
bindir=${project_dir}/bin

if [ $# -lt 3 ]
then
  echo "<Robject name> <non-transposed keytype (i.e. subject/probe)> <transposed keytype (i.e. subject/probe)>"
  exit 1
fi

robject=$1
notranspose=$2
withtranspose=$3

echo Converting ${robject}
sed "s/ROBJECT/${robject}/g" ${bindir}/import_rdata.template | sed "s/NOTRANSPOSE/${notranspose}/g" | sed "s/WITHTRANSPOSE/${withtranspose}/g" > import_rdata.r
R --no-save< import_rdata.r
sed "s/,/\t/" tmp_${notranspose}.csv > ${robject}_${notranspose}.tsv
sed "s/,/\t/" tmp_${withtranspose}.csv > ${robject}_${withtranspose}.tsv
exit 0
