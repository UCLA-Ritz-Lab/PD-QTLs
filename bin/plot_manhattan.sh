#!/bin/bash

if [ $# -lt  2 ]
then
  echo "Usage: <peg study (peg1|peg2)> <me_qtl_type (cis|trans)>" 
  exit 1
fi

peg=$1
me_qtl_type=$2

grep GSA- ${peg}_${me_qtl_type}_input.txt |cut -f3 > gsa.tmp

R --no-save << END
library(qqman)
gsa<-read.table("gsa.tmp",header=F)[,1]
gwasResults<-read.table("${peg}_${me_qtl_type}_input.txt",header=T,sep='\t')
cutoff<-gwasResults[10,4]
pdf("${peg}_${me_qtl_type}_output.pdf")
#png("${peg}_${me_qtl_type}_output.png",width=1024,height=768)
manhattan(gwasResults,highlight=gsa,annotatePval=cutoff,annotateTop=F,main="meQTLs for ${peg} of type ${me_qtl_type}")
dev.off()
END

