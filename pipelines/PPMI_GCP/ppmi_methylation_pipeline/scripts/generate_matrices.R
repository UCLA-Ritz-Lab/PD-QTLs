library(ChAMP)

csvfile=snakemake@input$sample_sheet
dir<-snakemake@input$idat
beta_col_file<-snakemake@output$beta_col
beta_row_file<-snakemake@output$beta_row
beta_val_file<-snakemake@output$beta_val
pd_col_file<-snakemake@output$pd_col
pd_row_file<-snakemake@output$pd_row
pd_val_file<-snakemake@output$pd_val

if (file.exists(csvfile)){
  print(csvfile)
  #print('exists')
  myImport<- champ.import(dir,arraytype="EPIC")
  write.table(colnames(myImport$beta),file=beta_col_file,row.names=F,col.names=F,quote=F,sep='\t')
  write.table(rownames(myImport$beta),file=beta_row_file,row.names=F,col.names=F,quote=F,sep='\t')
  write.table(myImport$beta,file=beta_val_file,row.names=F,col.names=F,quote=F,sep='\t')
  write.table(colnames(myImport$pd),file=pd_col_file,row.names=F,col.names=F,quote=F,sep='\t')
  write.table(rownames(myImport$pd),file=pd_row_file,row.names=F,col.names=F,quote=F,sep='\t')
  write.table(myImport$pd,file=,pd_val_file,row.names=F,col.names=F,quote=F,sep='\t')
}
