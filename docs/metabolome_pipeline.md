# Importing rawdata

## Obtain files from https://app.box.com/folder/170401803918?utm_source=trans

## Convert whole genome Rdata objects to TSV files

Note that the top level directory where this Git repo is downloaded will be denoted as [repo_root]. Under [repo_root], you should see directories such as docs, rawdata, and bin.

Under [repo_root]/rawdata/Metabolomics, run R:
	load('c18_keyVar_link_cc.RData')
	write.table(c18_keyVar_link_cc,file='covariates_c18.txt',sep='\t',quote=F,row.names=F,col.names=T)
	load('hilic_keyVar_link_cc.RData')
	write.table(hilic_keyVar_link_cc,file='covariates_hilic.txt',sep='\t',quote=F,row.names=F,col.names=T)
	load('PEG_hilic_quant_combat_PC_LOD_cc.t.RData')
	write.table(PEG_hilic_quant_combat_PC_cc.t,file='metabolites_hilic.txt',sep='\t',quote=F,row.names=F,col.names=T)
	load('PEG_c18_quant_combat_LOD_cc.t.RData')
	write.table(PEG_c18_quant_combat_cc.t,file='metabolites_c18.txt',sep='\t',quote=F,row.names=F,col.names=T)

	
Under [repo_root]/rawdata/Metabolomics:

```
../../bin/clip.py $'\t' 'Sample.ID,pegid' < covariates_c18.txt | sed '1d' > sample_map_c18_import.tsv
cut -f1 sample_map_c18_import.tsv > a
sed 's/\t/,/g' metabolites_c18.txt |sed '1d' > b
paste a b |sed s'/^/c18\t/' > c18_import.tsv
cut -f1-2,4-99 covariates_hilic.txt |../../bin/clip.py $'\t' 'Sample.ID,pegid' |  sed '1d' > sample_map_hilic_import.tsv
cut -f1 sample_map_hilic_import.tsv  > a
sed 's/\t/,/g' metabolites_hilic.txt |sed '1d' > b
paste a b |sed 's/^/hilic\t/'> hilic_import.tsv
../../bin/import_metabolome_data.sh

```
Under [repo_root]/rawdata/Metabolomics:
```
head -n1 metabolites_c18.txt  > a
head -n1 metabolites_hilic.txt  > b
paste a b > metabolite_list.txt
```
Under [repo_root]/rawdata/merge_metabolome/peg1:
```
../../../bin/fetch_raw_metabolome_matrices_peg.sh 1 peg1cases 1 | gzip -c - > raw_merge.txt.gz
### GKC CHECK HERE
cp ../../merge/peg1cases/snplist.txt .
zcat raw_merge.txt.gz |cut -f1|sed '1d' > subjectlist.txt
cp ../../Metabolomics/metabolite_list.txt .
gunzip -c raw_merge.txt.gz  | ../../../bin/process_merge.py snplist.txt cols subjectlist.txt rows metabolite_list.txt cols 6
../../../bin/transpose_float 437 8 < covariates.tsv > covariates_t.tsv
cat genotypes.csv |sed 's/\,/\t/g' | ../../../bin/transpose_float 437 263704 > genotypes_t.tsv
mv methylation.csv metabolites.csv
cat metabolites.csv |sed 's/\,/\t/g' | ../../../bin/transpose_float 437 5135 > metabolites_t.tsv
```

Under [repo_root]/rawdata/merge_metabolome/peg2:
```
../../../bin/fetch_raw_metabolome_matrices_peg.sh 2 209 | gzip -c - > raw_merge.txt.gz
cp ../../merge/peg1/snplist.txt .
zcat raw_merge.txt.gz |cut -f1|sed '1d' > subjectlist.txt
cp ../../Metabolomics/metabolite_list.txt .
gunzip -c raw_merge.txt.gz  | ../../../bin/process_merge.py snplist.txt cols subjectlist.txt rows metabolite_list.txt cols 6
../../../bin/transpose_float 206 8 < covariates.tsv > covariates_t.tsv
cat genotypes.csv |sed 's/\,/\t/g' | ../../../bin/transpose_float 206 263704 > genotypes_t.tsv
mv methylation.csv metabolites.csv
cat metabolites.csv |sed 's/\,/\t/g' | ../../../bin/transpose_float 206 5135 > metabolites_t.tsv
```
	
## Run matrix_eqtl

### all analysis

in [repo_root]/rawdata/merge_metabolome/peg[1|2]:

        R --no-save < ../../../bin/matrix_eqtl.r 1>cis_trans.out 2>cis_trans.err &

post process the results in [repo_root]/rawdata/merge_metabolome/peg[1|2]:
	sed 's/_/\t/' all_eqtls.raw | sed '1d' | sed 's/^/all\t/' > all_eqtls.txt


load the results in [repo_root]/rawdata/merge_metabolome:

        ../../bin/load_all_mqtls.sh

### Ontology enrichment analysis

For BED file in [repo_root]/rawdata/merge_metabolome/peg[1|2]:
        ../../../bin/get_eqtls.sh m_qtls all peg1 |sed '1d' |sed 's/^/chr/' > cis_eqtls.bed
        ../../../bin/get_eqtls.sh m_qtls all peg2 |sed '1d' |sed 's/^/chr/' > cis_eqtls.bed

### PD enrichment analysis

#### PEG1
in [repo_root]/rawdata/merge_metabolome/peg1:

#####all
        ../../../bin/fetch_hypergeometric_param.sh m_qtls all peg1

So among 263704 QC passed SNPs tested for association on PEG1, 201 of these were deemed as significant mQTLs.  Computing a p-value for the hypergeometic test where alternative hypothesis is observing 6 or more mQTLs among the 11920 PD SNPs.

	#all
 	> phyper(6-1,201,263704,11920,lower.tail=F)
	[1] 0.8945533

#### PEG2
in [repo_root]/rawdata/merge_metabolome/peg2:

#####all
        ../../../bin/fetch_hypergeometric_param.sh m_qtls all peg2

So among 272674 QC passed SNPs tested for association on PEG2, 148 of these were deemed as significant mQTLs.  Computing a p-value for the hypergeometic test where alternative hypothesis is observing 6 or more mQTLs among the 12120 PD SNPs.

	#all
 	> phyper(6-1,148,272674,12120,lower.tail=F)
	[1] 0.6469277

## getting hotspots

in [repo_root]/rawdata/merge_metabolome/peg[1|2]:
	
	../../../bin/get_hotspots.py < all_eqtls.raw |sed '1d' | sort -k2 -g -r | sed 's/^/metab\tpeg1\t/' > metab_hotspots.txt
	../../../bin/get_hotspots.py < all_eqtls.raw |sed '1d' | sort -k2 -g -r | sed 's/^/metab\tpeg2\t/' > metab_hotspots.txt


## getting meQTL metQTL overlap

in [repo_root]/results/eqtls
	../../bin/fetch_methyl_metab_qtl_overlap.sh peg1 > meqtl_metqtl_overlap_peg1.txt
	../../bin/fetch_methyl_metab_qtl_overlap.sh peg2 > meqtl_metqtl_overlap_peg2.txt
