dat<-read.table(file='methylation.csv',sep=',',header=T)
write.table(t(dat),file='methylation_t.tsv',quote=F,sep='\t',row.names=T,col.names=F)
