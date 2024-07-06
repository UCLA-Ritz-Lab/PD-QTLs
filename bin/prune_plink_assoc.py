#!/usr/bin/python3

import sys

def get_snp_set(filename):
	delimiter = '\t'
	snp_set = set()
	#header 
	f_handle = open(filename,'r')
	line = f_handle.readline()
	for line in f_handle:
		line = line.strip()
		tokens = line.split(delimiter)
		snp_set.add(tokens[1])
	f_handle.close()
	return snp_set

def print_snp_set(snp_set):
	for info in snp_set:
		print(info)

def main():
	if(len(sys.argv)<1):
		print("Usage: prune_plink_assoc.py [eqtl_file]")
		exit(1)
	eqtl_file = sys.argv[1]
	#plink_assoc_file = sys.argv[2]

	# some constants
	delimiter = '\t'

	snp_set = get_snp_set(eqtl_file)
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
	delimiter='\s+'
	tokens = line.split()
	#print(tokens[0]+','+tokens[1])
	for line in sys.stdin:
		line = line.rstrip()
		tokens = line.split()
		snp_id = tokens[1]
		if snp_id in snp_set:
			print(line)

#		header_tokens = []
#		f_covariates.write(subjectlist[linenum])
#		f_genotypes.write(subjectlist[linenum])
#		f_methylation.write(subjectlist[linenum])
#		if (linenum==0):
#			header_tokens = line.split(delimiter)
#			colnum = 0
#			for header_token in header_tokens:
#				if(colnum>=covariate_start_col and colnum<=covariate_end_col):
#					f_covariates.write(delimiter+header_token)
#				colnum+=1
#			f_covariates.write('\n')
#			for snp in snplist:
#				f_genotypes.write(','+snp)
#			f_genotypes.write('\n')
#			for probe in probelist:
#				f_methylation.write(','+probe)
#			f_methylation.write('\n')
#		else:
#			print('data line: '+str(linenum))
#			tokens = line.split(delimiter)
#			colnum = 0
#			for token in tokens:
#				if(colnum>=covariate_start_col and colnum<=covariate_end_col):
#					f_covariates.write(delimiter+token)
#					
#				elif(colnum==gwas_data_col):
#					f_genotypes.write(','+token)
#				elif(colnum==meth_data_col):
#					f_methylation.write(','+token)
#				colnum+=1
#			f_covariates.write('\n')
#			f_genotypes.write('\n')
#			f_methylation.write('\n')
#		linenum+=1

if __name__=="__main__":
	main()
