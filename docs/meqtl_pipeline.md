# Importing rawdata

## Clone a Git repo containing useful scripts and add to path

	git clone git@github.com:gchen98/data_tools.git

If repo is under $HOME directory:

	export PATH=$PATH:$HOME/data_tools


## Obtain files from https://app.box.com/folder/158643387254

## Convert whole genome methylation Rdata objects to TSV files

Note that the top level directory where this Git repo is downloaded will be denoted as [repo_root]. Under [repo_root], you should see directories such as docs, rawdata, and bin.

Download the files under the Data folder from the Box website.  Unzip files so that the contents are under [repo_root]/rawdata. 

On Hoffman cluster, under [repo_root]/rawdata/Methylation, for each of the RData files, you can run under R (after editing the script to either point to peg1 or peg2):

	source('champ_pipeline.r')

Files of interest:

	combat_peg1_by_subject.txt
	combat_peg1_by_cpg.txt
	peg1_dmp.txt
	imputed_peg2_by_subject.txt
	imputed_peg2_by_cpg
	gsea.txt

Copy to [repo_root]/rawdata/Methylation.txt

Prepare files for meta-analysis of DMPs

In R:
	setwd('[repo_root]/rawdata/DMP_meta')
	df_peg1<-read.table('../Methylation/peg1_dmp.txt',header=T,sep=' ',fill=TRUE)
	colnames(df_peg1)<-gsub('X1_to_0.','',colnames(df_peg1))
	df_peg1a<-cbind(rownames(df_peg1),df_peg1)
	colnames(df_peg1a)[1] = 'CpG'
	write.csv(df_peg1a,file='../DMP_meta/peg1.txt',row.names=F,quote=F)

In [repo_root]/rawdata/DMP_meta:

	../metal < metal.script

In R:

	df_dmp_meta<-read.table('../DMP_meta/DMP_OUT1.metal',sep='\t',header=T)
	


rewrite the header based on mapping files to prepare for Matrix EQTL:
	
	rewrite_header.py mapping_peg1.csv ',' combat_peg1_by_cpg.txt ' ' > ../MatrixEQTL/peg1/mapped_peg1_by_cpg.txt
	rewrite_header.py mapping_peg2.csv ',' imputed_peg2_by_cpg.txt  ' ' > ../MatrixEQTL/peg2/mapped_peg2_by_cpg.txt

## Extract critical covariates

Under [repo_root]/rawdata/Covariates:

	clip.py ',' Ethnicity,ExternalDNACode,Female,Age,Mono,Gran,CD4T,NK,CD8_naive,CD8pCD28nCD45RAn,PlasmaBlast,PDstudyDiseaseNumeric < all_covariates_peg1.csv  > covariates_peg1.csv
	clip.py ',' EthnicityOriginal,PEGID,Female,Age,PlasmaBlast,CD8pCD28nCD45RAn,CD4T,NK,Mono,Gran < all_covariates_peg2.csv | sed 's/$/,1/' > covariates_peg2.csv 

Edit the fields (ExternalDNACode) to (peg_id) in covariates_peg2.csv
Edit the fields (PEGID,EthnicityOriginal,1) to (peg_id,Ethnicity,PDstudyDiseaseNumeric) in covariates_peg2.csv
Reorder the columns:

	intersect_by_header.py ',' covariates_peg1.csv covariates_peg2.csv
	

Merge ancestry with clipped covariates:
	R --no-save<merge.r 
	cp peg1_*by_covariate.txt ../MatrixEQTL/peg1
	cp peg2_*by_covariate.txt ../MatrixEQTL/peg2

Copy covariates to Genetics folder:
	cp peg1_ancestry_by_subject.txt  ../Genetics/peg_covariates.txt
	sed '1d' peg2_ancestry_by_subject.txt  >> ../Genetics/peg_covariates.txt


##  Process whole genome VCF files 

Under [repo_root]/rawdata/Genetics:

Get a sample manifest:

	plink2 --vcf PEG_PD.phased.vcf.gz --make-pgen -out results/qc/qc1/PEG_PD_RAW

Generate pheno.txt from covariates:

	R --no-save< R/make_pheno.r

Filter on covariates

	plink2 --pfile results/qc/qc1/PEG_PD_RAW --pheno pheno.txt --make-pgen --require-pheno DISEASE -out results/qc/qc1/PEG_PD_COV_FILT


Compute HW stats:

	plink2 --pfile results/qc/qc1/PEG_PD_COV_FILT --hardy --snps-only --autosome --keep-founders --thin-indiv-count 100 --thin-count 100000 --out results/qc/qc1/PEG_PD_HW_ALL_POP

Prepare libraries in R:

	source("../R/install_packages.R")
	source("../R/loadplink2HW.R")
	source("../R/popstrat.R")
	source("../R/gwas1.R")

To get the HW Ternary plot in R:
	fname = "./results/qc/qc1/PEG_PD_HW_ALL_POP.hardy"
	data<-loadplink2Rmat(fname)
	HWTernaryPlot(data[,5:7],region=0,markercol = rgb(0,0,1,0.03),vertexlab = c("AA","AB","BB"))

and the HW QQ plot:

	format_hwdata2hwpack <- function(data) {
	  x <- data[, 5:7]  # Select relevant columns
	  colnames(x) <- c("AA", "AB", "BB")  # Set column names
	  
	  # Convert each column to numeric
	  numeric_matrix <- do.call(cbind, lapply(x, function(column) as.numeric(unlist(column))))
	  
	  return(numeric_matrix)  # Return the numeric matrix
	}
	data2 <- format_hwdata2hwpack(data)
	HWQqplot(data2,nsim=10,logplot=TRUE, main="Observed QQ HWE-Null")
See the allele frequency spectrum:

	plink2 --vcf PEG_PD.phased.vcf.gz \
	--freq alt1bins=0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.20,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.30,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.40,0.41,0.42,0.43,0.44,0.45,0.46,0.47,0.48,0.49,0.50,0.51,0.52,0.53,0.54,0.55,0.56,0.57,0.58,0.59,0.60,0.61,0.62,0.63,0.64,0.65,0.66,0.67,0.68,0.69,0.70,0.71,0.72,0.73,0.74,0.75,0.76,0.77,0.78,0.79,0.80,0.81,0.82,0.83,0.84,0.85,0.86,0.87,0.88,0.89,0.90,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,1.00 \
	--snps-only \
	--autosome \
	--keep-founders \
	--thin-count 100000 \
	--out results/qc/qc1/PEG_PD_freq_ALLPOP_keepall

In R:
	fname='results/qc/qc1/PEG_PD_freq_ALLPOP_keepal.afreq.alt1.bins'
	data = loadplink2Rmat(fname)
	barplot(unlist(data['OBS_CT']),names.arg=unlist(data['BIN_START']),ylim=c(0,1000),xlab="ALT1 FREQ",ylab="AUTOSOMAL SNP COUNT")

Remove SNPs with MAF <.01 call rates <.1 and out of HWE at p<1e-7

	plink2 --pfile results/qc/qc1/PEG_PD_COV_FILT --geno 0.1  --maf 0.01 -hwe 0.0000001 --export vcf --out results/qc/qc2/PEG_PD.phased.filtered

Generate Chr bgzip files under results/qc/qc2

	./loop_filter.sh

will generate new files in directory chunks

	#bgzip PEG_PD.phased.filtered.vcf
	#tabix -p vcf PEG_PD.phased.filtered.vcf.gz

Impute on imputation server:
	https://imputationserver.sph.umich.edu/#!run/imputationserver2

Save results in [repo_root]/rawdata/Genetics/imputed.  Enter that direc.

	../unzip_imputed.sh

Merge chr vcfs into a single one.

	ls chr*.dose.vcf.gz >filelist.txt
	bcftools concat -f filelist.txt -O z -o chr_all.dose.vcf.gz
	


#	Impute with BEAGLE:

#	java -jar beagle/beagle.27Feb25.75f.jar gt=results/qc/qc2/PEG_PD.phased.filtered.vcf out=results/qc/qc2/PEG_PD.phased.filtered.imputed

Convert back to plink2:

	plink2 --vcf imputed/chr_all.dose.vcf.gz dosage=DS --extract-if-info "R2 >= 0.8" --exclude-if-info "MAF < 0.01" --make-pgen --pheno pheno.txt -out results/qc/qc2/PEG_PD_IMPUTED
	#plink2 --vcf results/qc/qc2/PEG_PD.phased.filtered.imputed.vcf.gz --make-pgen --pheno pheno.txt -out results/qc/qc2/PEG_PD_IMPUTED

LD pruning:

	plink2 --pfile results/qc/qc2/PEG_PD_IMPUTED --rm-dup force-first --make-pgen --out results/qc/qc2/PEG_PD_dedup
	plink2 --pfile results/qc/qc2/PEG_PD_dedup --indep-pairwise 100 0.8 --out results/qc/qc3/PEG_PD_prune

Generate new dataset of independent variants:

	plink2 --pfile results/qc/qc2/PEG_PD_dedup --make-pgen --extract results/qc/qc3/PEG_PD_prune.prune.in --out results/qc/qc3/PEG_PD_LE


Generate kinship matrix using KING algorithm:

#	plink2 --pfile PEG_PD_LE --make-king square0 --maf 0.1 --out results/qc/qc3/PEG_PD_KING

Compute Principal Components:

	plink2 --pfile results/qc/qc3/PEG_PD_LE --freq counts --pca biallelic-var-wts --out results/qc/qc3/PEG_PD_pca
	sed 's/^#//' results/qc/qc3/PEG_PD_LE.psam  > psam.tmp
	sed 's/^#//' results/qc/qc3/PEG_PD_pca.eigenvec > evec.tmp

Generate PCA plots in R:

	setwd('[repo_root]/rawdata/Genetics')
	source('../R/popstrat.R')
	source('../R/popstrat_plot2d.r')
	popstrat_plot2d(df_merge,c('PC1','PC2'),Aprefix='ETHNICITY')
	outliers<-df_pc[df_pc$PC1<  (-0.2),]
	write.table(outliers$IID,file='results/qc/qc3/pc1_outliers.txt',quote=F,row.names=F,col.names=F)


Run second pass of Principal Components:

	plink2 --pfile results/qc/qc3/PEG_PD_LE --remove results/qc/qc3/pc1_outliers.txt --freq counts --pca biallelic-var-wts --out results/qc/qc3/PEG_PD_pca2
	sed 's/^#//' results/qc/qc3/PEG_PD_pca2.eigenvec > evec.tmp

Generate new version with samples filtered out.

	plink2 --pfile results/qc/qc2/PEG_PD_IMPUTED --remove results/qc/qc3/pc1_outliers.txt --make-pgen --out results/qc/qc3/PEG_PD_IMPUTED2

Chunk into chromosomes:

	for i in {1..22}; do
  	  plink2 --pfile results/qc/qc3/PEG_PD_IMPUTED2 --chr $i --export Av --out chunks/PEG_PD_CHR_$i
	done


Prepare for GWAS:

In R:

	df_covar<-read.table('covar.txt',header=T)
	df_eth<-read.table('evec.tmp',header=T,sep='\t')
	df_merge<-merge(df_covar,df_eth,by=c('IID'))
	write.table(df_merge,file='covar_evec.txt',row.names=F,quote=F)

Run GWAS:

	plink2 --pfile results/qc/qc3/PEG_PD_IMPUTED2 --glm omit-ref hide-covar no-firth --pheno-name DISEASE --covar covar_evec.txt --covar-name AGE,FEMALE,PC1,PC2 --vif 1000 --covar-variance-standardize -pfilter 1 --out results/gwas/PEG_PD
	sed 's/^#//' results/gwas/PEG_PD.DISEASE.glm.logistic > ../GWAS_meta/gwas_peg.txt

Prep for meta-analysis in R:
	setwd('[repo_root]/rawdata/Genetics/')
	df_peg<-read.table('GWAS_meta/gwas.txt',sep='\t',header=T)	
	df_nalls<-read.table('../Nalls_GWAS/PDNallsP005.csv',header=T,sep=',')
	df_nalls$ID<-paste(df_nalls$chr,df_nalls$position,sep=':')
	df_nalls$N<-df_nalls$ncase+df_nalls$ncontrol
	write.csv(df_nalls,file='../GWAS_meta/gwas_nalls.txt',row.names=F,quote=F)
	df_merged<-merge(df_peg,df_nalls,by.x=c('CHROM','POS'),by.y=c('chr','position'))

Run meta in [repo_root]/rawdata/GWAS_meta:
	../metal < metal.script

Sort:
	df_meta<-read.table('GWAS_OUT1.metal',sep='\t',header=T)
	write.table(ordered,file='sorted_meta.txt',quote=F,row.names=F)



Make a dataframe file from plink2:

	plink2 --vcf PEG_PD.phased.filtered.imputed.vcf.gz --export Av --out PEG_PD.dataframe

Post-process this dataframe to split into two files:

	./process_plink_dataframe.py < PEG_PD.dataframe.traw

rewrite the header based on mapping files to prepare for Matrix EQTL:

	rewrite_header.py mapping_peg.csv ',' PEG_PD.geno.txt ' ' > ../MatrixEQTL/mapped_peg_by_snp.txt

In MatrixEQTL make symlinks
In <repo_root>/rawdata/MatrixEQTL/peg1 for example:

	ln -s ../mapped_peg_by_snp.txt

## Intersect columns:

	intersect_by_header.py ' ' peg1_ancestry_by_covariate.txt mapped_peg_by_snp.txt mapped_peg1_by_cpg.txt

Filter on DMPs

	filter_on_whitelist.py ../dmp_cpgs.txt  ' ' < intersected_mapped_peg1_by_cpg.txt > intersected_mapped_peg1_by_cpg_dmp.txt

## Run MatrixEQTL

	


produces

	plink2.log  plink2.pgen  plink2.psam  plink2.pvar



## Import Ethnicity Information

In <repo_root>/rawdata/Ethnicity

    sed 's/"//g' GC_dataset.csv | ../../bin/clip.py ',' 'Pegid,Race_gwas,Race_based_on_Q,K1a,K2a,K3a,K4a' > gc_import.csv

Note Cynthia's note on two subjects with missing data.  The revised file is now gc_import_revised.csv

    sql_pd_qtl < ../../bin/create_ethnicity.sql



## Import PEG1 covariates into DB

Loading PEG1 covariates that include blood cell count variables

Generate a SQL import friendly file from the PEG1 CSV file by running. Clip only relevant variables with command:

	../../bin/clip.py ',' ExternalDNACode,SampleID,Ethnicity,Female,Age,RFvoteHispanic,PDstudyParkinsonsDisease,Mono,Gran,CD4T,NK,CD8.naive,CD8pCD28nCD45RAn,PlasmaBlast < PEG1_cov.csv > peg1_import.csv

Generate a draft SQL script for import by running:

	../../bin/make_create_table.py peg1_covariates ',' < peg1_import.csv  > ../../bin/create_peg1_covariates.sql

Edit the draft SQL script with vi:

	vi ../../bin/create_peg1_covariates.sql

Edit accordingly, changing SampleID (the sample ID on the microarray) to varchar(255) and adding a key for this so that we can resolve SampleIDs to PEG IDs (ExternalDNACode). Some fields should be changed to float, and PDstudyParkinsonsDisease can be a tinyint.  My version looks like:

```
use pd_qtl;
drop table if exists  peg1_covariates ;
create table  peg1_covariates  (peg_id  varchar(25), sample_id  varchar(25), Age  float, Ethnicity  text, PlasmaBlast  float, CD8pCD28nCD45RAn  float, CD8_naive  float, CD4T  float, NK  float, Mono  float, Gran  float, PDstudyParkinsonsDisease  tinyint, RFvoteHispanic  float,  primary key (peg_id), unique key index_sample_id(sample_id));
load data infile '/home/garyc/analysis/PD-QTLs/rawdata/Methylation/peg1_import.csv' into table  peg1_covariates  fields terminated by ',' ignore 1 lines;
```

## Import PEG2 covariates into DB

	sql_pd_qtl < ../../bin/fetch_peg2_link_sort.sql  | sed '1d' > GWAS.EWAS_209_link_sort.txt

	sql_pd_qtl < ../../bin/create_peg2_covariates.sql 

### Generate genotypes import file for MySQL

Before we do that let's generate the import file for the PEGID to GWASID mapping:

	../../bin/clip.py ',' 'Pegid,GWAS_ID,DUP' < Lill\ CRG\ GWAS\ link.csv > peg_id_gwas_id_mapping_import.csv

	../../bin/make_create_table.py peg_id_gwas_id_mapping ',' < peg_id_gwas_id_mapping_import.csv > ../../bin/create_peg_id_gwas_id_mapping.sql

Edit the SQL script. My version reads:

```
use pd_qtl;
drop table if exists  peg_id_gwas_id_mapping ;
create table  peg_id_gwas_id_mapping  (peg_id varchar(25), gwas_id varchar(25), dup tinyint,  key index_peg_id(peg_id), key index_gwas_id(gwas_id));
load data infile '/home/garyc/analysis/PD-QTLs/rawdata/Genetics/peg_id_gwas_id_mapping_import.csv' into table peg_id_gwas_id_mapping  fields terminated by ',' ignore 1 lines;
```

Import the file:

	sql_pd_qtl < ../../bin/create_peg_id_gwas_id_mapping.sql

In the directory [repo_root]/rawdata/Genetics, run:

	zcat PEG_PD.phased.vcf.gz |../../bin/parse_vcf.py 

This will create three files: snpinfo.txt, subjects.txt and genotypes.txt.  To import these three files into the DB proceed with the following three steps:



Generate SQL scripts for snpinfo.txt:

	../../bin/make_create_table.py snpinfo '\t' < snpinfo.txt  > ../../bin/create_snpinfo.sql

Edit to:

```
use pd_qtl;
drop table if exists  snpinfo ;
create table  snpinfo  (snpid varchar(255),chrom varchar(10),pos int unsigned, ref_allele varchar(10),alt_allele varchar(10),primary key(snpid),unique key index_position(chrom,pos));
load data infile '/home/garyc/analysis/PD-QTLs/rawdata/Genetics/snpinfo.txt' into table  snpinfo  fields terminated by '\t' ignore 0 lines;
```

Now import into DB:

	sql_pd_qtl < ../../bin/create_snpinfo.sql



Generate SQL scripts for subjects.txt:

	../../bin/make_create_table.py gwas_subjects '\t' < subjects.txt > ../../bin/make_gwas_subjects.sql

Edit make_gwas_subjects.sql to:

```
use pd_qtl;
drop table if exists  gwas_subjects ;
create table  gwas_subjects  (seq int auto_increment, gwas_id  varchar(255),  primary key ( seq ),unique key index_gwas_id(gwas_id));
load data infile '/home/garyc/analysis/PD-QTLs/rawdata/Genetics/subjects.txt' into table  gwas_subjects  fields terminated by '\t' ignore 0 lines(gwas_id);
```
 
Import:

	sql_pd_qtl < ../../bin/make_gwas_subjects.sql



# Use PLINK2 to filter GWAS dataset to keep SNPs that are 1) not in methylation probes, 2) >=5% genotyped 3) >=5% MAF 4) HWE p>=1e-7

##Generate a set of PLINK files based on filtering criteria:

### Prep subject filter files for plink

	echo "select gwas_id,gwas_id from peg_id_gwas_id_mapping as a, peg1_covariates as b where a.peg_id=b.peg_id and b.PDstudyParkinsonsDisease=1 order by gwas_id"|sql_pd_qtl |sed '1d' > GWAS.EWAS_peg1cases_link_sort.txt
	echo "select gwas_id,gwas_id from peg_id_gwas_id_mapping as a, peg1_covariates as b where a.peg_id=b.peg_id and b.PDstudyParkinsonsDisease=0 order by gwas_id"|sql_pd_qtl |sed '1d' > GWAS.EWAS_peg1controls_link_sort.txt
	mv GWAS.EWAS_209_link_sort.txt GWAS.EWAS_peg2cases_link_sort.txtt

### PEG1 cases

	plink2 --vcf PEG_PD.phased.vcf.gz --double-id --vcf-require-gt --geno 0.05  --maf 0.05 -hwe 0.0000001 --keep GWAS.EWAS_peg1cases_link_sort.txt --indiv-sort file GWAS.EWAS_peg1cases_link_sort.txt --make-bed --out PEG.phased.peg1cases

### PEG1 controls

	plink2 --vcf PEG_PD.phased.vcf.gz --double-id --vcf-require-gt --geno 0.05  --maf 0.05 -hwe 0.0000001 --keep GWAS.EWAS_peg1controls_link_sort.txt --indiv-sort file GWAS.EWAS_peg1controls_link_sort.txt --make-bed --out PEG.phased.peg1controls

### PEG2 cases

	plink2 --vcf PEG_PD.phased.vcf.gz --double-id --vcf-require-gt --geno 0.05  --maf 0.05 -hwe 0.0000001 --keep GWAS.EWAS_peg2cases_link_sort.txt --indiv-sort file GWAS.EWAS_peg2cases_link_sort.txt --make-bed --out PEG.phased.peg2cases

## Annotate the famfiles to have disease status to do GWAS

### PEG1

```
cp PEG.phased.580.fam PEG.phased.580.fam.bak
../../bin/annotate_famfile.sh PEG.phased.580.fam.bak |sed '1d' > PEG.phased.580.fam
../../bin/fetch_gwas_covar.sh PEG.phased.580.fam.bak |sed '1d' |sed 's/NULL/-9/g' > PEG.phased.580.cov
```


## Run GWAS and load
```
#plink2 --bfile PEG.phased.580 --glm --fam PEG.phased.580.fam --covar PEG.phased.580.cov
#./../bin/load_gwas_results.sh plink2.PHENO1.glm.logistic.hybrid
plink --bfile PEG.phased.580 --logistic --allow-no-sex  --covar PEG.phased.580.cov
#TO DO
#./../bin/load_gwas_results.sh plink2.PHENO1.glm.logistic.hybrid

```

## Clump SNPs

```
plink --bfile PEG.phased.580 --clump plink.assoc.logistic  --clump-p1 0.001 --clump-p2 0.01 --clump-r2 0.5 --clump-kb 250
```

## Run shell scripts to filter out probes near SNPs. From [repo_root]/rawdata/Methylation:

	sql_pd_qtl < ../../bin/fetch_meth_exclusion_windows.sql |sed '1d' > METH_exclude_sql.txt

##Generate a new set of PLINK files based on CpG range exclusion for SNPs

### PEG1 cases
	plink2 --bfile PEG.phased.peg1cases --exclude 'range' METH_exclude_sql.txt --make-bed --out PEG.phased.peg1cases.methex
### PEG1 controls
	plink2 --bfile PEG.phased.peg1controls --exclude 'range' METH_exclude_sql.txt --make-bed --out PEG.phased.peg1controls.methex
### PEG2 cases
	plink2 --bfile PEG.phased.peg2cases --exclude 'range' METH_exclude_sql.txt --make-bed --out PEG.phased.peg2cases.methex

# Generate subject major dataset with additive genotypes

## Using PLINK to recode to additive dosage:

	plink2 --bfile PEG.phased.580.methex --recode A --out PEG.phased.580.methex.AD

or

### PEG1 cases
	plink2 --bfile PEG.phased.peg1cases --recode A --out PEG.phased.peg1cases.AD
### PEG1 controls
	plink2 --bfile PEG.phased.peg1controls --recode A --out PEG.phased.peg1controls.AD
### PEG2
	plink2 --bfile PEG.phased.peg2cases --recode A --out PEG.phased.peg2cases.AD


##Make SNP map file

	plink2 --bfile PEG.phased.580.methex --recode bimbam --out PEG.phased.580.methex.bimbam

or

### PEG1 cases
	plink2 --bfile PEG.phased.peg1cases --recode bimbam --out PEG.phased.peg1cases.bimbam
### PEG1 controls
	plink2 --bfile PEG.phased.peg1controls --recode bimbam --out PEG.phased.peg1controls.bimbam
### PEG2 cases
	plink2 --bfile PEG.phased.peg2cases --recode bimbam --out PEG.phased.peg2cases.bimbam


## Concatenate subject's genotypes into a single CSV genotype string:

	cut -f1-6 --complement PEG.phased.580.methex.AD.raw |sed 's/\t/,/g' > b
	cut -f2  PEG.phased.580.methex.AD.raw |paste - b > gwas_methex_additive_genotypes.txt

or

### PEG1 cases
	cut -f1-6 --complement PEG.phased.peg1cases.AD.raw |sed 's/\t/,/g' > b
	cut -f2  PEG.phased.peg1cases.AD.raw |paste - b > gwas_additive_genotypes_peg1cases.txt
### PEG1 controls
	cut -f1-6 --complement PEG.phased.peg1controls.AD.raw |sed 's/\t/,/g' > b
	cut -f2  PEG.phased.peg1controls.AD.raw |paste - b > gwas_additive_genotypes_peg1controls.txt
### PEG2 cases
	cut -f1-6 --complement PEG.phased.peg2cases.AD.raw |sed 's/\t/,/g' > b
	cut -f2  PEG.phased.peg2cases.AD.raw |paste - b > gwas_additive_genotypes_peg2cases.txt

run the script in [repo_root]/bin:

### PEG1 cases
	./create_gwas_subject_genotypes.sh peg1cases
### PEG1 controls
	./create_gwas_subject_genotypes.sh peg1controls
### PEG2 cases
	./create_gwas_subject_genotypes.sh peg2cases


## Dump a merge of PEG1 Methylation data

To get a dataset for all subjects that have covariates, genotypes, and methylation data, we run a join on MySQL. 


### PEG1 cases
In [repo_root]/rawdata/merge/peg1cases:
	../../../bin/fetch_raw_matrices_peg.sh 1 peg1cases 1 | gzip -c - > raw_merge.txt.gz

### PEG1 controls
In [repo_root]/rawdata/merge/peg1controls:
	../../../bin/fetch_raw_matrices_peg.sh 1 peg1controls 0 | gzip -c - > raw_merge.txt.gz

### PEG2 cases
In [repo_root]/rawdata/merge/peg2cases:
	../../../bin/fetch_raw_matrices_peg.sh 2 peg2cases 1 | gzip -c - > raw_merge.txt.gz

## Extract the SNP lists

Next we must provide the raw_merge.txt.gz file with some metadata including a SNP list, subject list, and probe list.

in [repo_root]/rawdata/Genetics:

### PEG1 cases
	head -n1 PEG.phased.peg1cases.AD.raw |cut -f1-6 --complement > ../merge/peg1cases/snplist.txt
### PEG1 controls
	head -n1 PEG.phased.peg1controls.AD.raw |cut -f1-6 --complement > ../merge/peg1controls/snplist.txt
### PEG2 cases
	head -n1 PEG.phased.peg2cases.AD.raw |cut -f1-6 --complement > ../merge/peg2cases/snplist.txt


## Get a list of the subjects:

### PEG1 cases
in [repo_root]/rawdata/merge/peg1cases:
	gunzip -c raw_merge.txt.gz |cut -f1|sed '1d' > subjectlist.txt
### PEG1 controls
in [repo_root]/rawdata/merge/peg1controls:
	gunzip -c raw_merge.txt.gz |cut -f1|sed '1d' > subjectlist.txt
### PEG2 cases
in [repo_root]/rawdata/merge/peg2cases:
	gunzip -c raw_merge.txt.gz |cut -f1|sed '1d' > subjectlist.txt


##Get a list of the probes 

in [repo_root]/rawdata/Methylation:

### PEG1 cases
	cut -f1 datMethPEG1t_probe.tsv  > ../merge/peg1cases/probelist.txt
### PEG1 controls
	cut -f1 datMethPEG1t_probe.tsv  > ../merge/peg1controls/probelist.txt
### PEG2 cases
	cut -f1 datMethPEG2t_probe.tsv  > ../merge/peg2cases/probelist.txt

## Run the script to process the merge 

in [repo_root]/rawdata/merge/[peg1cases|peg1controls|peg2cases]:
	gunzip -c raw_merge.txt.gz  | ../../../bin/process_merge.py snplist.txt cols subjectlist.txt rows probelist.txt rows 13

## Transpose the three matrices

### PEG1 cases
in [repo_root]/rawdata/merge/peg1cases:
```
cat covariates.tsv |../../../bin/transpose_float 319 11 > covariates_t.tsv
cat genotypes.csv |sed 's/\,/\t/g' | ../../../bin/transpose_float 319 267770 > genotypes_t.tsv 
cat methylation.csv |sed 's/\,/\t/g' | ../../../bin/transpose_float 319 485512 > methylation_t.tsv 
```
### PEG1 controls
in [repo_root]/rawdata/merge/peg1controls:
```
cat covariates.tsv |../../../bin/transpose_float 232 11 > covariates_t.tsv
cat genotypes.csv |sed 's/\,/\t/g' | ../../../bin/transpose_float 232 267770 > genotypes_t.tsv 
cat methylation.csv |sed 's/\,/\t/g' | ../../../bin/transpose_float 232 485512 > methylation_t.tsv 
```
### PEG2 cases
in [repo_root]/rawdata/merge/peg2cases:
```
cat covariates.tsv |../../../bin/transpose_float 209 11 > covariates_t.tsv
cat genotypes.csv |sed 's/\,/\t/g' | ../../../bin/transpose_float 209 272674 > genotypes_t.tsv 
cat methylation.csv |sed 's/\,/\t/g' | ../../../bin/transpose_float 209 485512 > methylation_t.tsv 

```


## cis trans analysis

### Generate a list of probes that do not contain SNPs

These are eligible for eQTL analysis for cis trans. To run

in [repo_root]/rawdata/Methylation:
	sql_pd_qtl < ../../bin/fetch_nosnp_probes.sql |sed '1d' > nosnp_probes.txt

### To get a methylation dataset that excludes probes with SNPs in it:

in [repo_root]/rawdata/merge/[peg1cases|peg1controls|peg2cases]:
	../../../bin/filter_probes.py whitelist ../../Methylation/nosnp_probes.txt < methylation_t.tsv > methylation_nosnp_probes_t.tsv

### To get a gene map and snp map file 

in [repo_root]/rawdata/merge:
```
sql_pd_qtl < ../../bin/fetch_nosnp_probes_and_map.sql > gene_map.txt
sql_pd_qtl < ../../bin/fetch_gwas_map.sql|sed 's/\t/_/' >snp_map.txt
```

### Run the analysis 


in [repo_root]/rawdata/merge/[peg1cases|peg1controls|peg2cases]:

	R --no-save < ../../../bin/matrix_eqtl.r 1>cis_trans.out 2>cis_trans.err &

post process the results:

```
#sed 's/_/\t/' cis_eqtls.raw | sed '1d' | sed 's/^/cis\t/' > cis_eqtls.txt
#sed 's/_/\t/' trans_eqtls.raw | sed '1d' | sed 's/^/trans\t/' > trans_eqtls.txt
sed 's/_/\t/' all_eqtls.raw | sed '1d' | sed 's/^/all\t/' > all_eqtls.txt
cat cis_eqtls.raw | sed '1d' | sed 's/^/cis\t/' > cis_eqtls.txt
cat trans_eqtls.raw | sed '1d' | sed 's/^/trans\t/' > trans_eqtls.txt

```

load the results:

	../../../bin/load_all_meqtls.sh

Dump all results:

In <repo_root>/results/eqtls:
	../../bin/fetch_all_meqtls_distinct.sh > cases_only_cistrans_meqtl.txt
	head -n1 cases_only_cistrans_meqtl.txt > allcases.header
	grep ^cis cases_only_cistrans_meqtl.txt  | sed '1d'|sort -k14 -g | cat allcases.header - > cases_only_cis_meqtls.txt 


Stratify results by cohort and cis/trans
	head -n1 cases_only_cistrans_meqtl.txt | cut -f2,5-6,8-11 > peg1cases.header
	head -n1 cases_only_cistrans_meqtl.txt | cut -f2,5-6,13-16 > peg2cases.header

	grep ^cis cases_only_cistrans_meqtl.txt |grep peg1cases | cut -f2,5-6,8-11 > tmp; cat peg1cases.header tmp >  peg1cases_cis.txt
	grep ^cis cases_only_cistrans_meqtl.txt |grep peg2cases |  cut -f2,5-6,13-16 >tmp; cat peg2cases.header tmp > peg2cases_cis.txt
	grep ^trans cases_only_cistrans_meqtl.txt |grep peg1cases | cut -f2,5-6,8-11 > tmp; cat peg1cases.header tmp >peg1cases_trans.txt
	grep ^trans cases_only_cistrans_meqtl.txt |grep peg2cases | cut -f2,5-6,13-16 > tmp; cat peg2cases.header tmp >peg2cases_trans.txt


### Running METAL for meta analysis

In <repo_root>/results/eqtls:
	metal < metal.script
	cat cis1.metal|sed '1d'|sed 's/^/cis\t/' |sed 's/,/\t/' |cut -f1,2,3,7,8 > metal.import 
	cat trans1.metal|sed '1d'|sed 's/^/trans\t/' |sed 's/,/\t/' |cut -f1,2,3,7,8 >> metal.import 
	../../bin/load_metal_results.sh

### PD enrichment analysis

In [repo_root]/rawdata/Genetics:
	../../bin/clip.py $'\t' "SNP,Beta_all_studies,SE_all_studies,P_all_studies" < PDNallsP005_TableS2.csv > pdnalls_table_s2_import.tsv


load pd probes in [repo_root]/bin:

	sql_pd_qtl < create_pd_probes.sql

load PD snps in [repo_root]/bin:

	sql_pd_qtl < create_pd_snps.sql

#### Overlap with existing GWAS SNPs

In <repo_root>/results/eqtls:
	../../bin/fetch_all_meqtls_pd_annot.sh > pdsnps_meqtls_peg.txt

#### Hypergeometric test
```
> phyper(15, 150, 400-150, 50)
[1] 0.1549789
> fisher.test(matrix(c(15, 50-15, 150-15, 400 - 50 - 150 + 15), nr = 2),
              alternative = "less")$p.value
[1] 0.1549789

https://stats.stackexchange.com/questions/288081/use-fishers-exact-test-or-a-hypergeometric-test
```

in [repo_root]/rawdata/merge/[peg1cases|peg2cases|peg1controls]:

	../../../bin/row2col.sh < snplist.txt > snplist_t.txt

#####cis
	../../../bin/fetch_hypergeometric_param.sh me_qtls cis peg1cases

in R:
```
gwas_snps<- 267770
gwas_eqtls<- 22968
pd_snps<- 11920
pd_qtls<- 1467
phyper(pd_qtls-1,gwas_eqtls,gwas_snps-gwas_eqtls,pd_snps,lower.tail=F)
[1] 2.687785e-45
#odds ratio
((pd_qtls/pd_snps)/(1-pd_qtls/pd_snps))/((gwas_eqtls/gwas_snps)/(1-gwas_eqtls/gwas_snps))
fisher.test(matrix(c(pd_qtls-1, pd_snps-pd_qtls-1, gwas_eqtls-pd_qtls-1, gwas_snps-gwas_eqtls-pd_snps+pd_qtls-1), nr = 2), alternative = "greater")$p.value

> phyper(15, 150, 400-150, 50)
[1] 0.1549789
> fisher.test(matrix(c(15, 50-15, 150-15, 400 - 50 - 150 + 15), nr = 2),
              alternative = "less")$p.value
```

## OLD VERSION
So among 263704 QC passed SNPs tested for association on PEG1, 33430 of these were deemed as significant cis meQTLs.  Computing a p-value for the hypergeometic test where alternative hypothesis is observing 1966 or more cis meQTLs among the 11920 PD SNPs.

	> phyper(1966-1, 33430, 263704-33430, 11920,lower.tail=F)
	[1] 2.681329e-35


#####trans
	../../../bin/fetch_hypergeometric_param.sh me_qtls trans peg1

```
gwas_snps<-267770
gwas_eqtls<- 28746
pd_snps<- 11920
pd_qtls<- 1059
phyper(pd_qtls-1,gwas_eqtls,gwas_snps-gwas_eqtls,pd_snps,lower.tail=F)
1
```

## OLD VERSION
So among 263704 QC passed SNPs tested for association on PEG1, 27904 of these were deemed as significant trans meQTLs.  Computing a p-value for the hypergeometic test where alternative hypothesis is observing 1161 or more trans meQTLs among the 11920 PD SNPs.

	> phyper(1161-1, 27904, 263704-27904, 11920,lower.tail=F)
	[1] 0.9990464
	

#### PEG2

in [repo_root]/rawdata/merge/peg2:
#####cis
	../../../bin/fetch_hypergeometric_param.sh me_qtls cis peg2

```
gwas_snps<- 272674
gwas_eqtls<- 18982
pd_snps<- 12120
pd_qtls<- 1212
phyper(pd_qtls-1,gwas_eqtls,gwas_snps-gwas_eqtls,pd_snps,lower.tail=F)
[1] 2.843806e-37
#odds ratio
((pd_qtls/pd_snps)/(1-pd_qtls/pd_snps))/((gwas_eqtls/gwas_snps)/(1-gwas_eqtls/gwas_snps))

```


### OLD VERSION
So among 272674 QC passed SNPs tested for association on PEG2, 19037 of these were deemed as significant cis meQTLs.  Computing a p-value for the hypergeometic test where alternative hypothesis is observing 1208 or more cis meQTLs among the 12120 PD SNPs.

	#cis
	> phyper(1208-1, 19037, 272674-19037, 12120,lower.tail=F)
	[1] 5.014762e-36

#####trans
```
gwas_snps<- 272674
gwas_eqtls<- 29679
pd_snps<- 12120
pd_qtls<- 1233
phyper(pd_qtls-1,gwas_eqtls,gwas_snps-gwas_eqtls,pd_snps,lower.tail=F)
[1] 0.9954382
```
### OLD VERSION
	../../../bin/fetch_hypergeometric_param.sh me_qtls trans peg2
So among 272674 QC passed SNPs tested for association on PEG2, 29999 of these were deemed as significant trans meQTLs.  Computing a p-value for the hypergeometic test where alternative hypothesis is observing 1293 or more trans meQTLs among the 12120 PD SNPs.

	#trans
	> phyper(1293,29999,272674-29999,12120,lower.tail=F)
	[1] 0.8823499

### Overlap with BIOS repo on meQTLs (molgenis)

In <repo_root>/rawdata/external_qtl:

	zcat 2015_09_02_cis_meQTLsFDR0.05-CpGLevel.txt.gz | ../../bin/clip.py $'\t' PValue,SNPName,ProbeName,CisTrans,AlleleAssessed,OverallZScore,FDR > molgenis_import.txt
	sql_pd_qtl < ../..bin/create_molgenis.sql

In <repo_root>/results/eqtls:

	../../bin/fetch_all_meqtls_molgenis_annot.sh > molgenis_meqtls_peg.txt


### GKC

### Ontology enrichment analysis

In [repo_root]/results/great:
	#../../bin/get_eqtl_positions.sh me_qtls cis peg1cases |sed '1d' | sed 's/^/chr/' > peg1cases_cis_eqtls.bed
	#../../bin/get_eqtl_positions.sh me_qtls cis peg2cases |sed '1d' | sed 's/^/chr/' > peg2cases_cis_eqtls.bed
	#../../bin/get_eqtl_positions.sh me_qtls cis all|sed '1d' |sed 's/^/chr/' > cis_eqtls.bed
	awk  '{print $3"\t"$4"\t"$4}' ../eqtls/cases_only_cistrans_meqtl.txt |sed '1d' |sed 's/^/chr/' > all_cis_eqtls.bed

### Pruning SNPs

Make sure we install TwoSampleMR package from https://github.com/MRCIEU/TwoSampleMR. We can now create a pruned snplist for running the clump_data function in peg[1|2]:

	cut -f1-2 cis_eqtls.bed |sed 's/\t/:/' > ldlinkr_snplist.txt

Run the R script to generate pruned SNPs for peg[1|2]:

	R --no-save < ../../../bin/prune_snplist.r

[GO term enrichments PEG1](http://www.caseyandgary.com:8099/~garyc/pd_qtl/go_enrichment_peg1.html)
[GO term enrichments PEG2](http://www.caseyandgary.com:8099/~garyc/pd_qtl/go_enrichment_peg2.html)

### trans hotspots

In <repo_root>/results/eqtls:
	#../../bin/get_hotspots.py < overlap_cistrans_meqtls_peg1-sorted.txt |sort -k2 -g -r  > overlap_trans_meqtls_peg1-sorted.txt
	../../bin/get_hotspots.py < cases_only_cistrans_meqtl.txt |sort -k4 -g  > cases_only_trans_meqtl.txt


# PEG2 pipeline

Repeat PEG1 with peg1 string replaced by peg2

# Manhattan plots

In <repo_root>/results/manhattan, generate the input files for the plots:

	../../bin/fetch_manhattan_plot_input.sh peg1 cis > peg1_cis_input.txt

Do the same for peg2 and trans as well.


# coloc analysis

In <repo_root>/rawdata/bin:
 
For cis_trans analysis, set cis p-value threshold in matrix_eqtl.r to one to output all results.

Rerun matrix_eqtl.r

Rename cis_eqtls.txt to cis_eqtls_all_p.txt

	LC_ALL=C sort -b  -k2 cis_eqtls_all_p.txt > cis_eqtls_all_p_sorted.txt.sorted

In <repo_root>/rawdata/Genetics:

	LC_ALL=C sort -b -k3 plink2.PHENO1.glm.logistic.hybrid > plink2.PHENO1.glm.logistic.hybrid.sorted

Join the two files

	LC_ALL=C join -1 3 -2 2 plink2.PHENO1.glm.logistic.hybrid.sorted ../merge/peg1/cis_eqtls_all_p_sorted.txt|  grep ADD  > coloc_merged.txt
	cat coloc_merged.txt |sort -t\  -k2n -k17d -k3n > coloc_merged_sorted.txt

# Running ewas and gwas screen for features correlated to PD status in rawdata/merge/peg1all

	../../../bin/run_screen.py disease.txt [genotypes_t.tsv|methylation_t.tsv] [gwas_pvalues.txt|ewas_pvalues.txt]
o
# selecting top 20 SNPs and CPGs from the analysis above from <repo_root>/results/metqtls/

	./fetch_top_results.sh
