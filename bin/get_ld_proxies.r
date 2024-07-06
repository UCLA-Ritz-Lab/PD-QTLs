library(LDlinkR)
token<-'e4b3b190d410'
dataset<-read.table('snp_positions.txt',header=T,sep='\t')
results<-LDproxy_batch(dataset$SNP,pop="EUR",r2d="r2",token=token,append = TRUE)
