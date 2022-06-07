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

### PEG1

	plink2 --vcf PEG_PD.phased.vcf.gz --double-id --vcf-require-gt --geno 0.05  --maf 0.05 -hwe 0.0000001 --keep GWAS.EWAS_580_link_sort.txt --indiv-sort file GWAS.EWAS_580_link_sort.txt --make-bed --out PEG.phased.580

### PEG2

	plink2 --vcf PEG_PD.phased.vcf.gz --double-id --vcf-require-gt --geno 0.05  --maf 0.05 -hwe 0.0000001 --keep GWAS.EWAS_209_link_sort.txt --indiv-sort file GWAS.EWAS_209_link_sort.txt --make-bed --out PEG.phased.209

## Run shell scripts to filter out probes near SNPs. From [repo_root]/rawdata/Methylation:

	sql_pd_qtl < ../../bin/fetch_meth_exclusion_windows.sql |sed '1d' > METH_exclude_sql.txt

##Generate a new set of PLINK files based on CpG range exclusion for SNPs

### PEG1
	plink2 --bfile PEG.phased.580 --exclude 'range' METH_exclude_sql.txt --make-bed --out PEG.phased.580.methex
### PEG2
	plink2 --bfile PEG.phased.209 --exclude 'range' METH_exclude_sql.txt --make-bed --out PEG.phased.209.methex

# Generate subject major dataset with additive genotypes

## Using PLINK to recode to additive dosage:

	plink2 --bfile PEG.phased.580.methex --recode A --out PEG.phased.580.methex.AD

or

### PEG1
	plink2 --bfile PEG.phased.580 --recode A --out PEG.phased.580.AD
### PEG2
	plink2 --bfile PEG.phased.209 --recode A --out PEG.phased.209.AD


##Make SNP map file

	plink2 --bfile PEG.phased.580.methex --recode bimbam --out PEG.phased.580.methex.bimbam

or

### PEG1
	plink2 --bfile PEG.phased.580 --recode bimbam --out PEG.phased.580.bimbam
### PEG2
	plink2 --bfile PEG.phased.209 --recode bimbam --out PEG.phased.209.bimbam


## Concatenate subject's genotypes into a single CSV genotype string:

	cut -f1-6 --complement PEG.phased.580.methex.AD.raw |sed 's/\t/,/g' > b
	cut -f2  PEG.phased.580.methex.AD.raw |paste - b > gwas_methex_additive_genotypes.txt

or

### PEG1
	cut -f1-6 --complement PEG.phased.580.AD.raw |sed 's/\t/,/g' > b
	cut -f2  PEG.phased.580.AD.raw |paste - b > gwas_additive_genotypes_580.txt
### PEG2
	cut -f1-6 --complement PEG.phased.209.AD.raw |sed 's/\t/,/g' > b
	cut -f2  PEG.phased.209.AD.raw |paste - b > gwas_additive_genotypes_209.txt

run the script in [repo_root]/bin:

### PEG1
	./create_gwas_subject_genotypes.sh 580
### PEG2
	./create_gwas_subject_genotypes.sh 209


## Dump a merge of PEG1 Methylation data

To get a dataset for all subjects that have covariates, genotypes, and methylation data, we run a join on MySQL. 

### PEG1
In [repo_root]/rawdata/merge/peg1:
	../../../bin/fetch_raw_matrices_peg.sh 1 580 | gzip -c - > raw_merge.txt.gz
### PEG2
In [repo_root]/rawdata/merge/peg2:
	../../../bin/fetch_raw_matrices_peg.sh 2 209 | gzip -c - > raw_merge.txt.gz

Next we must provide the raw_merge.txt.gz file with some metadata including a SNP list, subject list, and probe list.

Get a list of the SNPs in [repo_root]/rawdata/Genetics:

	head -n1 PEG.phased.580.methex.AD.raw |cut -f1-6 --complement > ../merge/peg1/snplist.txt

or 


### PEG1
	head -n1 PEG.phased.580.AD.raw |cut -f1-6 --complement > ../merge/peg1/snplist.txt
### PEG2
	head -n1 PEG.phased.209.AD.raw |cut -f1-6 --complement > ../merge/peg2/snplist.txt


Get a list of the subjects:

### PEG1
in [repo_root]/rawdata/merge/peg1:
	gunzip -c raw_merge.txt.gz |cut -f1|sed '1d' > subjectlist.txt
### PEG2
in [repo_root]/rawdata/merge/peg2:
	gunzip -c raw_merge.txt.gz |cut -f1|sed '1d' > subjectlist.txt

st
Get a list of the probes in [repo_root]/rawdata/Methylation:

### PEG1
	cut -f1 datMethPEG1t_probe.tsv  > ../merge/peg1/probelist.txt
### PEG2
	cut -f1 datMethPEG2t_probe.tsv  > ../merge/peg2/probelist.txt

Run the script to process the merge 

### PEG1
in [repo_root]/rawdata/merge/peg1:
	gunzip -c raw_merge.txt.gz  | ../../../bin/process_merge.py snplist.txt cols subjectlist.txt rows probelist.txt rows
### PEG2
in [repo_root]/rawdata/merge/peg2:
	gunzip -c raw_merge.txt.gz  | ../../../bin/process_merge.py snplist.txt cols subjectlist.txt rows probelist.txt rows

Transpose the three matrices
### PEG1
in [repo_root]/rawdata/merge/peg1:
```
cat covariates.tsv |../../../bin/transpose_float 551 11 > covariates_t.tsv
#cat genotypes.csv |sed 's/\,/\t/g' | ../../../bin/transpose_float 551 64326 > genotypes_t.tsv 
cat genotypes.csv |sed 's/\,/\t/g' | ../../../bin/transpose_float 551 263704 > genotypes_t.tsv 
cat methylation.csv |sed 's/\,/\t/g' | ../../../bin/transpose_float 551 485512 > methylation_t.tsv 
```
### PEG1
in [repo_root]/rawdata/merge/peg2:
```
cat covariates.tsv |../../../bin/transpose_float 209 11 > covariates_t.tsv
#cat genotypes.csv |sed 's/\,/\t/g' | ../../../bin/transpose_float 209 64326 > genotypes_t.tsv 
cat genotypes.csv |sed 's/\,/\t/g' | ../../../bin/transpose_float 209 272674 > genotypes_t.tsv 
cat methylation.csv |sed 's/\,/\t/g' | ../../../bin/transpose_float 209 485512 > methylation_t.tsv 
```

## cis trans analysis

### Generate a list of probes that do not contain SNPs

These are eligible for eQTL analysis for cis trans. To run

in [repo_root]/rawdata/Methylation:
	sql_pd_qtl < ../../bin/fetch_nosnp_probes.sql |sed '1d' > nosnp_probes.txt

To get a methylation dataset that excludes probes with SNPs in it:

#### PEG1
in [repo_root]/rawdata/merge/peg1:
	../../../bin/filter_probes.py whitelist ../../Methylation/nosnp_probes.txt < methylation_t.tsv > methylation_nosnp_probes_t.tsv
#### PEG2
in [repo_root]/rawdata/merge/peg2:
	../../../bin/filter_probes.py whitelist ../../Methylation/nosnp_probes.txt < methylation_t.tsv > methylation_nosnp_probes_t.tsv

### To get a gene map and snp map file 

in [repo_root]/rawdata/merge:
```
sql_pd_qtl < ../../bin/fetch_nosnp_probes_and_map.sql > gene_map.txt
sql_pd_qtl < ../../bin/fetch_gwas_map.sql|sed 's/\t/_/' >snp_map.txt
```

### Run the analysis 

in [repo_root]/rawdata/merge/peg[1|2]:

```
R --no-save < ../../../bin/matrix_eqtl.r
sed 's/_/\t/' cis_eqtls.raw | sed '1d' | sed 's/^/cis\t/' > cis_eqtls.txt
sed 's/_/\t/' trans_eqtls.raw | sed '1d' | sed 's/^/trans\t/' > trans_eqtls.txt
sed 's/_/\t/' all_eqtls.raw | sed '1d' | sed 's/^/all\t/' > all_eqtls.txt
```

Load the results and get the BED file:

#### PEG1
	../../../bin/get_cis_eqtls.sh 1 |sed '1d' |sed 's/^/chr/' > cis_eqtls.bed
#### PEG2
	../../../bin/get_cis_eqtls.sh 2 |sed '1d' |sed 's/^/chr/' > cis_eqtls.bed

[GO term enrichments](http://www.caseyandgary.com:8099/~garyc/pd_qtl/go_enrichment.html)

[Full results from GREAT enrichment for cis-trans](http://www.caseyandgary.com:8099/~garyc/pd_qtl/cis_eqtls.bed)

### PD disease SNPs for enrichment analysis

Load in the SNPs 

#### PEG1
in [repo_root]/rawdata/merge/peg1:
	../../../bin/fetch_hypergeometric_param.sh 1
#### PEG2
in [repo_root]/rawdata/merge/peg2:
	../../../bin/fetch_hypergeometric_param.sh 2

So among 64,326 QC passed SNPs tested for association on PEG, 2,216 of these were deemed as significant meQTLs.  Computing a p-value for the hypergeometic test where alternative hypothesis is observing 120 or more meQTLs among the 2602 PD SNPs.

```
#cis
> phyper(120-1, 2216, 64326-2216, 2602,lower.tail=F)
[1] 0.000834084
> phyper(353-1,4170,64326-4170,2602,lower.tail=F)
[1] 2.172357e-40
#trans
> phyper(244-1,6605,64326-6605,2602,lower.tail=F)
[1] 0.9421406
> phyper(842-1,13596,64326-13596,2602,lower.tail=F)
[1] 2.272072e-42
#all
> phyper(157-1,3681,64326-3681,2602,lower.tail=F)
[1] 0.2540671
> phyper(948-1,8409,64326-8409,2602,lower.tail=F)
[1] 3.100249e-211
```
 

# PEG2 pipeline

Repeat PEG1 with peg1 string replaced by peg2


# Results

## Open issue

Covariates?

```
garyc@lupine:~/analysis/PD-QTLs/rawdata/merge$ head -n2 peg?/covariates.tsv
==> peg1/covariates.tsv <==
subject	Female	Age	RFvoteHispanic	PDstudyParkinsonsDisease	Mono	Gran	CD4T	NK	CD8_naive	CD8pCD28nCD45RAn	PlasmaBlast
10002AP40	0	61	0.153439	1	0.0883461	0.702761	0.108779	0.0183139	105.07	14.5551	2.07082

==> peg2/covariates.tsv <==
subject	Female	Age	RFvoteHispanic	PDstudyParkinsonsDisease	Mono	Gran	CD4T	NK	CD8_naive	CD8pCD28nCD45RAn	PlasmaBlast
80159WF40	0	71.2	NULL	1	0.048608	0.555363	0.0541631	0.146495	166.338	14.6923	1.79475
garyc@lupine:~/analysis/PD-QTLs/rawdata/merge$ `
```

## PEG2 replication

### 73026 of 109178 eQTLS from PEG1 are not replicated in PEG2
### 0 of 35908 eQTLS from PEG2 are not replicated in PEG1
```
garyc@lupine:~/analysis/PD-QTLs/results/eqtls$ grep -c NULL *txt
distinct_meqtls_peg1.txt:73026
distinct_meqtls_peg2.txt:0
overlap_meqtls_peg1.txt:0
overlap_meqtls_peg2.txt:0
garyc@lupine:~/analysis/PD-QTLs/results/eqtls$ wc -l distinct_meqtls_peg1.txt 
109178 distinct_meqtls_peg1.txt
garyc@lupine:~/analysis/PD-QTLs/results/eqtls$ wc -l distinct_meqtls_peg2.txt 
35908 distinct_meqtls_peg2.txt
garyc@lupine:~/analysis/PD-QTLs/results/eqtls$ 
```

## Top results from PEG1 or PEG2

```
garyc@lupine:~/analysis/PD-QTLs/results/eqtls$ head distinct_meqtls_peg*
==> distinct_meqtls_peg1.txt <==
me_qtl_type	snp_id	allele	gene	peg	pvalue	FDR	beta	peg	pvalue	FDR	beta
cis	GSA-rs1040961	G	cg17707870	peg1	3.19657356478992e-240	3.02266263743662e-232	0.463461	peg2	1.75221926908519e-73	9.04349479281351e-67	0.476207
cis	rs10010994	C	cg17858192	peg1	1.07333906461107e-224	5.07471800999251e-217	0.38231	peg2	2.87543180201319e-69	8.83254921602953e-63	0.380936
cis	rs10184015	A	cg02502145	peg1	2.35478400502351e-214	7.42223664073641e-207	0.458667	peg2	3.27733181336013e-71	1.46082830374007e-64	0.477257
cis	rs2532925	G	cg04145681	peg1	9.42822327913758e-211	2.22881920794558e-203	0.471904	peg2	2.27317819557558e-74	1.31125151237358e-67	0.458799
cis	exm2267473	G	cg09084244	peg1	2.65120419044779e-207	5.01392862592143e-200	0.422279	peg2	2.66697691072591e-68	5.23059017069958e-62	0.407247
cis	GSA-rs1035142	T	cg07227024	peg1	3.20562587043044e-203	5.05203741431135e-196	0.428486	peg2	5.38180254456733e-71	2.1989631919515e-64	0.434122
cis	GSA-rs10750097	G	cg12556569	peg1	1.21710492165548e-201	1.64412545315535e-194	0.385482	peg2	2.47609121473652e-70	9.33887867099282e-64	0.394837
cis	rs1043793	A	cg02502145	peg1	2.89813128784325e-199	3.42557154739125e-192	0.452316	peg2	7.18611377732221e-44	1.61625187488907e-38	0.429988
cis	rs8106375	G	cg22996768	peg1	7.34618381232944e-199	7.7183462184694e-192	0.472055	peg2	1.70179415323157e-55	9.11921341179653e-50	0.466453

==> distinct_meqtls_peg2.txt <==
me_qtl_type	snp_id	allele	gene	peg	pvalue	FDR	beta	peg	pvalue	FDR	beta
cis	GSA-rs4796640	G	cg25929399	peg1	7.80749012719841e-181	3.35578197621924e-174	0.438301	peg2	1.90294376752955e-111	1.86606770493654e-103	0.499874
cis	GSA-rs2294942	G	cg05704942	peg1	1.22762395772679e-183	7.25521600440376e-177	0.403242	peg2	7.78514468065094e-100	3.81714039970854e-92	0.421001
cis	rs1939015	G	cg10306192	peg1	1.90285543502298e-190	1.63575435080116e-183	0.492246	peg2	1.16352650586597e-91	3.803263963063e-84	0.530363
cis	rs2883456	C	cg11144103	peg1	5.06679537113219e-187	3.42223874350835e-180	0.429745	peg2	3.53248196898722e-91	8.66007003606234e-84	0.460901
cis	rs72660967	T	cg06961873	peg1	6.02441299455294e-156	1.03575495915106e-149	0.404299	peg2	4.4902824016553e-87	8.80653555690919e-80	0.437046
cis	rs3762352	C	cg24088508	peg1	2.59899202418852e-169	6.14398192883974e-163	-0.3777	peg2	7.15474810507206e-87	1.16935005329171e-79	-0.390806
cis	rs6142884	A	cg00704664	peg1	1.08394110417597e-153	1.68027677565248e-147	0.366204	peg2	2.50136142729234e-85	3.5041270999555e-78	0.388578
cis	GSA-rs35850196	A	cg06961873	peg1	2.71839833826451e-154	4.28417122490655e-148	0.402319	peg2	6.74026000425308e-85	8.26205543447358e-78	0.434157
cis	rs61137192	G	cg22337626	peg1	2.37168470250421e-141	1.54665668976338e-135	0.384497	peg2	5.95694195188303e-84	6.49056196583822e-77	0.40965
garyc@lupine:~/analysis/PD-QTLs/results/eqtls$ 
```
## Top overlapping results

```
garyc@lupine:~/analysis/PD-QTLs/results/eqtls$ head overlap_meqtls_peg*
==> overlap_meqtls_peg1.txt <==
me_qtl_type	snp_id	allele	gene	peg	pvalue	FDR	beta	peg	pvalue	FDR	beta
cis	GSA-rs1040961	G	cg17707870	peg1	3.19657356478992e-240	3.02266263743662e-232	0.463461	peg2	1.75221926908519e-73	9.04349479281351e-67	0.476207
cis	rs10010994	C	cg17858192	peg1	1.07333906461107e-224	5.07471800999251e-217	0.38231	peg2	2.87543180201319e-69	8.83254921602953e-63	0.380936
cis	rs10184015	A	cg02502145	peg1	2.35478400502351e-214	7.42223664073641e-207	0.458667	peg2	3.27733181336013e-71	1.46082830374007e-64	0.477257
cis	rs2532925	G	cg04145681	peg1	9.42822327913758e-211	2.22881920794558e-203	0.471904	peg2	2.27317819557558e-74	1.31125151237358e-67	0.458799
cis	exm2267473	G	cg09084244	peg1	2.65120419044779e-207	5.01392862592143e-200	0.422279	peg2	2.66697691072591e-68	5.23059017069958e-62	0.407247
cis	GSA-rs1035142	T	cg07227024	peg1	3.20562587043044e-203	5.05203741431135e-196	0.428486	peg2	5.38180254456733e-71	2.1989631919515e-64	0.434122
cis	GSA-rs10750097	G	cg12556569	peg1	1.21710492165548e-201	1.64412545315535e-194	0.385482	peg2	2.47609121473652e-70	9.33887867099282e-64	0.394837
cis	rs1043793	A	cg02502145	peg1	2.89813128784325e-199	3.42557154739125e-192	0.452316	peg2	7.18611377732221e-44	1.61625187488907e-38	0.429988
cis	rs8106375	G	cg22996768	peg1	7.34618381232944e-199	7.7183462184694e-192	0.472055	peg2	1.70179415323157e-55	9.11921341179653e-50	0.466453

==> overlap_meqtls_peg2.txt <==
me_qtl_type	snp_id	allele	gene	peg	pvalue	FDR	beta	peg	pvalue	FDR	beta
cis	GSA-rs4796640	G	cg25929399	peg1	7.80749012719841e-181	3.35578197621924e-174	0.438301	peg2	1.90294376752955e-111	1.86606770493654e-103	0.499874
cis	GSA-rs2294942	G	cg05704942	peg1	1.22762395772679e-183	7.25521600440376e-177	0.403242	peg2	7.78514468065094e-100	3.81714039970854e-92	0.421001
cis	rs1939015	G	cg10306192	peg1	1.90285543502298e-190	1.63575435080116e-183	0.492246	peg2	1.16352650586597e-91	3.803263963063e-84	0.530363
cis	rs2883456	C	cg11144103	peg1	5.06679537113219e-187	3.42223874350835e-180	0.429745	peg2	3.53248196898722e-91	8.66007003606234e-84	0.460901
cis	rs72660967	T	cg06961873	peg1	6.02441299455294e-156	1.03575495915106e-149	0.404299	peg2	4.4902824016553e-87	8.80653555690919e-80	0.437046
cis	rs6142884	A	cg00704664	peg1	1.08394110417597e-153	1.68027677565248e-147	0.366204	peg2	2.50136142729234e-85	3.5041270999555e-78	0.388578
cis	GSA-rs35850196	A	cg06961873	peg1	2.71839833826451e-154	4.28417122490655e-148	0.402319	peg2	6.74026000425308e-85	8.26205543447358e-78	0.434157
cis	rs61137192	G	cg22337626	peg1	2.37168470250421e-141	1.54665668976338e-135	0.384497	peg2	5.95694195188303e-84	6.49056196583822e-77	0.40965
cis	rs1865574	A	cg18584561	peg1	9.14738499362999e-160	1.66340724445189e-153	0.325022	peg2	7.81305755263808e-80	7.66165276376831e-73	0.351741
garyc@lupine:~/analysis/PD-QTLs/results/eqtls$
```

