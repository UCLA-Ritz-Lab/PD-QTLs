#!/bin/bash

project_dir='/home/garyc/analysis/PD-QTLs'
bindir=${project_dir}/bin

if [ $# -lt 3 ]
then
  echo "<infile> <delimiter> <tsv outfile>"
  exit 1
fi

infile=$1
delimiter=$2
outfile=$3

echo Converting ${infile} with delimiter ${delimiter} to output TSV ${outfile}
sed "s/INFILE/${infile}/g" ${bindir}/transpose_rmatrix.template | sed "s/DELIMITER/${delimiter}/g" | sed "s/OUTFILE/${outfile}/g"  > transpose_rmatrix.r
R --no-save< transpose_rmatrix.r
exit 0
