#!/usr/bin/python3

import sys

def main():
	fdr_threshold = 1
	delimiter = '\t'
	comma = ','
	header = False
	cpg_count = {}
	cpg_meta_z = {}
	cpg_meta_p = {}
	cpg_unique = {}
	#f_snpinfo = open('snpinfo.txt','w')
	#f_genotypes = open('genotypes.txt','w')
	for line in sys.stdin:
		line = line.rstrip()
		tokens = line.split(delimiter)
		if (len(tokens)>2 and tokens[1] != 'snp_id'):
			qtl_type = tokens[0]
			snp_gene = tokens[1]
			tokens2 = snp_gene.split(comma)
			snp_id = tokens2[0]
			gene = tokens2[1]
			meta_z = tokens[14]
			meta_p = tokens[15]
			snp=snp_id
			#snp=snp_id+"_"+allele
			#print("snp "+snp+" "+meta_z+" "+meta_p)
			cpg = tokens[3]
			if (qtl_type == 'trans'):
				if (snp not in cpg_unique):
					cpg_unique[snp]=set()
				cpg_set = cpg_unique[snp]
				cpg_set.add(cpg)
				cpg_unique[snp] = cpg_set
				cpg_count[snp] = len(cpg_set)
				cpg_meta_z[snp] = meta_z
				cpg_meta_p[snp] = meta_p
	print("snp\tprobe_count\tmeta_z-score\tmeta_p-value\ttrans_genes")
	for snp in cpg_count.keys():
		cpg_set = cpg_unique[snp]
		print(snp+"\t"+str(cpg_count[snp])+"\t"+str(cpg_meta_z[snp])+"\t"+str(cpg_meta_p[snp])+"\t"+str(cpg_set))
		

if __name__=="__main__":
	main()
