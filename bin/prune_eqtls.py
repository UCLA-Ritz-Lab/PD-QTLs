#!/usr/bin/python3

import sys

def get_snp_set(filename):
	#delimiter = '\t'
	snp_set = set()
	#header 
	f_handle = open(filename,'r')
	line = f_handle.readline()
	line.strip()
	#print(line)
	for line in f_handle:
		line = line.strip()
		#print(line)
		tokens = line.split()
		#print(str(len(tokens)))
		if(len(tokens)>2):
			snp_set.add(tokens[2])
	f_handle.close()
	return snp_set

def print_snp_set(snp_set):
	for info in snp_set:
		print(info)

def main():
	if(len(sys.argv)<1):
		print("Usage: prune_eqtls.py [plink_clump_file]")
		exit(1)
	plink_clump_file = sys.argv[1]
	#plink_clump_file = sys.argv[2]

	# some constants
	delimiter = '\t'

	snp_set = get_snp_set(plink_clump_file)
	#probelist = get_snp_set(probe_info_file,probe_info_orient)
	#print('snps')
	#print_snp_set(snp_set)
	#print('subjects')
	#print_snp_set(subjectlist)
	#print('probes')
	#print_snp_set(probelist)
	linenum = 0
	#f_covariates = open('covariates.tsv','w')
	#f_genotypes = open('genotypes.csv','w')
	#f_methylation = open('methylation.csv','w')
	line = sys.stdin.readline()
	line = line.rstrip()
	print(line)
	delimiter='\t'
	tokens = line.split(delimiter)
	#print(tokens[0]+','+tokens[1])
	for line in sys.stdin:
		line = line.rstrip()
		tokens = line.split(delimiter)
		snp_id = tokens[1]
		if snp_id in snp_set:
			print(line)


if __name__=="__main__":
	main()
