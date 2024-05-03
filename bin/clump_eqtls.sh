#!/bin/bash

if [ $# -lt 1 ]
then
  echo "Usage: <eqtl_file>"
  exit 1
fi

eqtl_file=$1

bin_dir='../../bin'
plink_dir='../../rawdata/Genetics'

${bin_dir}/prune_plink_assoc.py $eqtl_file < ${plink_dir}/plink.assoc.logistic > plink.assoc.logistic.pruned 

plink --bfile ${plink_dir}/PEG.phased.580 --clump plink.assoc.logistic.pruned  --clump-p1 0.0001 --clump-r2 0.001 --clump-kb 10000 1>out 2>err

${bin_dir}/prune_eqtls.py plink.clumped < $eqtl_file
