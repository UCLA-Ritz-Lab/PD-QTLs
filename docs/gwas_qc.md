# Notes for QC of raw genotype PLINK2 file

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
	


#	Impute with BEAGLE (if run locally):

	#java -jar beagle/beagle.27Feb25.75f.jar gt=results/qc/qc2/PEG_PD.phased.filtered.vcf out=results/qc/qc2/PEG_PD.phased.filtered.imputed

Convert back to plink2:

	plink2 --vcf imputed/chr_all.dose.vcf.gz dosage=DS --extract-if-info "R2 >= 0.8" --exclude-if-info "MAF < 0.01" --make-pgen --pheno pheno.txt -out results/qc/qc2/PEG_PD_IMPUTED
	#plink2 --vcf results/qc/qc2/PEG_PD.phased.filtered.imputed.vcf.gz --make-pgen --pheno pheno.txt -out results/qc/qc2/PEG_PD_IMPUTED

LD pruning:

	plink2 --pfile results/qc/qc2/PEG_PD_IMPUTED --rm-dup force-first --make-pgen --out results/qc/qc2/PEG_PD_dedup
	plink2 --pfile results/qc/qc2/PEG_PD_dedup --indep-pairwise 100 0.8 --out results/qc/qc3/PEG_PD_prune

Generate new dataset of independent variants:

	plink2 --pfile results/qc/qc2/PEG_PD_dedup --make-pgen --extract results/qc/qc3/PEG_PD_prune.prune.in --out results/qc/qc3/PEG_PD_LE


Generate kinship matrix using KING algorithm:

	#plink2 --pfile PEG_PD_LE --make-king square0 --maf 0.1 --out results/qc/qc3/PEG_PD_KING

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

