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

Generate SQL scripts for genotypes.txt:

	../../bin/make_create_table.py genotypes '\t' < genotypes.txt  > ../../bin/make_genotypes.sql

Edit make_genotypes.sql to:

```
use pd_qtl;
drop table if exists  genotypes ;
create table  genotypes  (snpid varchar(255), genotype_string mediumtext, primary key(snpid));
load data infile '/home/garyc/analysis/PD-QTLs/rawdata/Genetics/genotypes.txt' into table  genotypes  fields terminated by '\t' ignore 0 lines;
```

Import:

	sql_pd_qtl < ../../bin/make_genotypes.sql


# Use PLINK2 to filter GWAS dataset to keep SNPs that are 1) not in methylation probes, 2) >=5% genotyped 3) >=5% MAF 4) HWE p>=1e-7

##Generate a set of PLINK files based on filtering criteria:

	plink2 --vcf PEG_PD.phased.vcf.gz --double-id --vcf-require-gt --geno 0.05  --maf 0.05 -hwe 0.0000001 --keep GWAS.EWAS_580_link_sort.txt --indiv-sort file GWAS.EWAS_580_link_sort.txt --make-bed --out PEG.phased.580

## Run shell scripts to filter out probes near SNPs. From [repo_root]/rawdata/Methylation:

	sql_pd_qtl < ../../bin/fetch_meth_exclusion_windows.sql |sed '1d' > METH_exclude_sql.txt

##Generate a new set of PLINK files based on CpG range exclusion for SNPs

	plink2 --bfile PEG.phased.580 --exclude 'range' METH_exclude_sql.txt --make-bed --out PEG.phased.580.methex

# Generate subject major dataset with additive genotypes

## Using PLINK to recode to additive dosage:

	plink2 --bfile PEG.phased.580.methex --recode A --out PEG.phased.580.methex.AD

##Make SNP map file

 plink2 --bfile PEG.phased.580.methex --recode bimbam --out PEG.phased.580.methex.bimbam

## Concatenate subject's genotypes into a single CSV genotype string:

	cut -f1-6 --complement PEG.phased.580.methex.AD.raw |sed 's/\t/,/g' > b
cut -f2  PEG.phased.580.methex.AD.raw |paste - b > gwas_additive_genotypes.txt

## Create a SQL create table and import script in [repo_root]/bin:

	./make_create_table.py gwas_subject_genotypes '\t' < ../rawdata/Genetics/gwas_additive_genotypes.txt > create_gwas_subject_genotypes.sql

## Edit the SQL script. My version is:

```
use pd_qtl;
drop table if exists  gwas_subject_genotypes ;
create table  gwas_subject_genotypes  (gwas_id varchar(25), genotype_string mediumtext, primary key(gwas_id));
load data infile '/home/garyc/analysis/PD-QTLs/rawdata/Genetics/gwas_additive_genotypes.txt' into table  gwas_subject_genotypes  fields terminated by '\t' ignore 1 lines;
```

# Dump a merge of PEG1 Methylation data

To get a dataset for all subjects that have covariates, genotypes, and methylation data, we run a join on MySQL. In [repo_root]/rawdata/merge:

	sql_pd_qtl < ../../bin/fetch_raw_matrices.sql | gzip -c - > raw_merge.txt.gz

Next we must provide the raw_merge.txt.gz file with some metadata including a SNP list, subject list, and probe list.

Get a list of the SNPs in [repo_root]/rawdata/Genetics:

	head -n1 PEG.phased.580.methex.AD.raw |cut -f1-6 --complement > ../merge/snplist.txt


Get a list of the subjects in [repo_root]/rawdata/merge:

	gunzip -c raw_merge.txt.gz |cut -f1|sed '1d' > subjectlist.txt

Get a list of the probes in [repo_root]/rawdata/Methylation:

	cut -f1 datMethPEG1t_probe.tsv  > ../merge/probelist.txt

Run the script to process the merge in [repo_root]/rawdata/merge:

	gunzip -c raw_merge.txt.gz  | ../../bin/process_merge.py snplist.txt cols subjectlist.txt rows probelist.txt rows

Transpose the three matrices in [repo_root]/rawdata/merge:

```
cat covariates.tsv |../../bin/transpose_float 551 11 > covariates_t.tsv
cat genotypes.csv |sed 's/\,/\t/g' | ../../bin/transpose_float 551 64326 > genotypes_t.tsv 
cat methylation.csv |sed 's/\,/\t/g' | ../../bin/transpose_float 551 485512 > methylation_t.tsv 
```

# cis trans analysis

## Generate a list of probes that do not contain SNPs

These are eligible for eQTL analysis for cis trans. To run

	sql_pd_qtl < ../../bin/fetch_nosnp_probes.sql |sed '1d' > nosnp_probes.txt
To get a methylation dataset that excludes probes with SNPs in it:

	../../bin/filter_probes.py whitelist ../Methylation/nosnp_probes.txt < methylation_t.tsv > methylation_nosnp_probes_t.tsv

## To get a gene map and snp map file 

in [repo_root]/rawdata/merge:

```
sql_pd_qtl < ../../bin/fetch_nosnp_probes_and_map.sql > gene_map.txt
sql_pd_qtl < ../../bin/fetch_gwas_map.sql|sed 's/\t/_/' >snp_map.txt
```

## Run the analysis 

in [repo_root]/rawdata/merge:

```
R --no-save < ../../bin/matrix_eqtls.r
sed 's/_/\t/' cis_eqtls.raw | sed '1d' > cis_eqtls.txt
```

Load the results and get the BED file:

	sql_pd_qtl < ../../bin/get_cis_eqtls.sql |sed '1d' |sed 's/^/chr/' > cis_eqtls.bed

[GO term enrichments](http://www.caseyandgary.com:8099/~garyc/pd_qtl/go_enrichment.html)

[Full results from GREAT enrichment for cis-trans](http://www.caseyandgary.com:8099/~garyc/pd_qtl/cis_eqtls.bed)


