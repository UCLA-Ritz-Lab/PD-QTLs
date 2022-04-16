METH <-read.csv("HumanMethylation450_15017482_v.1.2.csv")

myvars<- c('IlmnID','CHR','MAPINFO')
METH_sub <- METH[myvars]

#remove missing mapinfo (based on genome, build 37), mainifest includes ~1000 rows at bottom for rs's
METH_sub <- METH_sub[complete.cases(METH_sub),]

#create start/end location, 20kb
METH_sub$start=METH_sub$MAPINFO-20000
METH_sub$end=METH_sub$MAPINFO+20000

#re-set those with - value to 1
METH_sub$start<-ifelse(METH_sub$start<0, 1, METH_sub$start)

#remove x/y chr (no genotype)
METH_sub<-METH_sub[!(METH_sub$CHR=='Y'|METH_sub$CHR=='X'),]

#order variables for plink (chr, start, end, id (cpg))
Meth2=METH_sub[,c(2,4,5,1)]

options(scipen=10) #remove scientific notation
write.table(Meth2, file="METH_exclude.txt", quote=FALSE, sep="\t", col.names=F, row.names=F)
options(scipen=0)
