# Importing rawdata

## Obtain files from https://app.box.com/folder/158643387254

## Convert whole genome methylation Rdata objects to TSV files

Note that the top level directory where this Git repo is downloaded will be denoted as [repo_root]. Under [repo_root], you should see directories such as docs, rawdata, and bin.

Download the files under the Data folder from the Box website.  Unzip files so that the contents are under [repo_root]/rawdata. In [repo_root]/bin, edit the files rdata2tsv.sh and tsv2sql.sh so that the project_dir variable points to the [repo_root].

Under [repo_root]/rawdata/Methylation, for each of the RData files, you can run:

	../../bin/rdata2tsv.sh

which will give you the arguments necessary for converting the RData objects into text files that the SQL database can parse on import.  Note that the first argument for this shell script is simply the RData file name without the .RObject extension.  The second argument is the key type for the rows of the RObject and the third argument is the key type for the columns.  So if rows corresponded to subjects and columns as probes, you can specify 'subject' for the second argument and 'probe' for the third argument.  Once you have converted each of the RObjects, you will see the same number of new files as original RObject files in the same directory.  You can then proceed towards importing these into SQL tables.

## Automatically importing TSV files into database (DB)

For any given TSV file, under [repo_root]/rawdata/Methylation, run 

	../../bin/tsv2sql.sh

which will provide necessary arguments for this shell script.  The first argument is simply the name of the new tsv file that was generated in the previous step.  The second column will be the index that represents the rows in the matrix in the tsv file.  So the second argument could be subject_id or probe_id for example.

To loop over all the tsv files in the current directory, a convenience script is provided in [repo_root]/bin/ called load_all_tsv.sh.  This script can just be simply run in the same folder as the tsv files are located in, without any arguments passed in.  

Under [repo_root]/rawdata/Methylation, simply run 

	../../bin/load_all_tsv.sh

## Import methylation annotation into DB

The Python script clip.py is useful for clipping only the necessary columns you want to import into the DB.  It also converts empty strings into MySQL friendly NULL strings (e.g. Null fields will be printed as \N).  For more info, under [repo_root]/rawdata/Methylation, run
 
	../../bin/clip.py

which will list arguments required for running the script. You can specify the delimiter and the list of columns. The actual input will be piped in via STDIN.


Under [repo_root]/rawdata/Methylation, run 

	../../bin/clip.py ',' IlmnID,CHR,MAPINFO,Probe_SNPs,Probe_SNPs_10 < HumanMethylation450_15017482_v.1.2.csv > illumina_annotation_import.csv

which would save the new file as illumina_annotation_import.csv in the current folder.

There is a convenience Python script that generates a SQL import script, which in turn needs to be slightly edited so that certain db table columns can have the correct data type. For more details, run 

	../../bin/make_create_table.py

to display help on arguments.  In this case run

	../../bin/make_create_table.py illumina_annotation ',' < illumina_annotation_import.csv  > ../../bin/create_illumina_annotation.sql

which will save a SQL script in the bin folder. Open this script with vi to edit the file accordingly :

	vi ../../bin/create_illumina_annotation.sql 

Certain things that should be edited: Columns such as numerics can be converted from text for float for example.  Also, make sure the load data infile statement points to the illumina_annotation_import.csv file. We can create indices too. Make sure CHR is changed to chrom, and its datatype is varchar(10).  MAPINFO should be a an int type.  The string ", unique key index_map(chrom,MAPINFO)" should be added right after the string "primary key ( IlmnID ),".  

My version looks like:

```
use pd_qtl;
drop table if exists  illumina_annotation ;
create table  illumina_annotation  (IlmnID  varchar(255), chrom varchar(25), mapinfo int unsigned,  primary key ( IlmnID ), unique key index_position(chrom,mapinfo));
load data infile '/home/garyc/analysis/PD-QTLs/rawdata/Methylation/illumina_annotation_import.csv' into table  illumina_annotation  fields terminated by ',' ignore 1 lines;
```

Now run the script as:

	sql_pd_qtl < ../../bin/create_illumina_annotation.sql

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

in [repo_root]/rawdata/merge/[peg1cases|peg2cases|peg1controls]:

	../../../bin/row2col.sh < snplist.txt > snplist_t.txt

#####cis
	../../../bin/fetch_hypergeometric_param.sh me_qtls cis peg1

So among 263704 QC passed SNPs tested for association on PEG1, 33430 of these were deemed as significant cis meQTLs.  Computing a p-value for the hypergeometic test where alternative hypothesis is observing 1966 or more cis meQTLs among the 11920 PD SNPs.

	> phyper(1966-1, 33430, 263704-33430, 11920,lower.tail=F)
	[1] 2.681329e-35

## NEW VERSION
> phyper(1467-1,22968,267770-22968,11920,lower.tail=F)
[1] 2.687785e-45

#####trans
	../../../bin/fetch_hypergeometric_param.sh me_qtls trans peg1

So among 263704 QC passed SNPs tested for association on PEG1, 27904 of these were deemed as significant trans meQTLs.  Computing a p-value for the hypergeometic test where alternative hypothesis is observing 1161 or more trans meQTLs among the 11920 PD SNPs.

	> phyper(1161-1, 27904, 263704-27904, 11920,lower.tail=F)
	[1] 0.9990464
	

#### PEG2

in [repo_root]/rawdata/merge/peg2:
#####cis
	../../../bin/fetch_hypergeometric_param.sh me_qtls cis peg2

So among 272674 QC passed SNPs tested for association on PEG2, 19037 of these were deemed as significant cis meQTLs.  Computing a p-value for the hypergeometic test where alternative hypothesis is observing 1208 or more cis meQTLs among the 12120 PD SNPs.

	#cis
	> phyper(1208-1, 19037, 272674-19037, 12120,lower.tail=F)
	[1] 5.014762e-36

#####trans
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
	../../bin/get_hotspots.py < cases_only_cistrans_meqtl.txt |sort -k2 -g -r > cases_only_trans_meqtl.txt


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





