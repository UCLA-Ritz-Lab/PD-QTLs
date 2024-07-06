#library(LDlinkR)
library(TwoSampleMR)
snplist<-read.table('ldlinkr_snplist.txt')
token<-'e4b3b190d410'
snplist_new<-clump_data(snplist,clump_kb=10000,clump_r2=0.001,clump_p1=1e-4,pop="EUR")
#LDproxy_batch(snplist,pop="EUR",r2d="r2",token=token,append = TRUE)
