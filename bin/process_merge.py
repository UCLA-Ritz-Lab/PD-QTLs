#!/usr/bin/python3

import sys

def get_info_list(filename,orient):
	delimiter = '\t'
	info_list = []
	f_handle = open(filename,'r')
	if (orient == 'rows'):
		for line in f_handle:
			line = line.strip()
			info_list.append(line)
	elif (orient == 'cols'):
		line = f_handle.readline()
		tokens = line.split(delimiter)
		for token in tokens:
			token = token.strip()
			info_list.append(token)
	f_handle.close()
	return info_list

def print_info_list(info_list):
	for info in info_list:
		print(info)

def main():
	# some constants
	delimiter = '\t'
	covariate_start_col = 2
	covariate_end_col = 12
	gwas_data_col = 14
	meth_data_col = 15
	if(len(sys.argv)<7):
		print("Usage: process_merge.py [snp_info file] [cols|rows] [subject_info_file] [cols|rows] [probe_info file] [cols|rows]")
		exit(1)
	snp_info_file = sys.argv[1]
	snp_info_orient = sys.argv[2]
	subject_info_file = sys.argv[3]
	subject_info_orient = sys.argv[4]
	probe_info_file = sys.argv[5]
	probe_info_orient = sys.argv[6]
	snplist = get_info_list(snp_info_file,snp_info_orient)
	subjectlist = get_info_list(subject_info_file,subject_info_orient)
	subjectlist.insert(0,'subject')
	probelist = get_info_list(probe_info_file,probe_info_orient)
	print('snps')
	#print_info_list(snplist)
	print('subjects')
	#print_info_list(subjectlist)
	print('probes')
	#print_info_list(probelist)
	linenum = 0
	f_covariates = open('covariates.tsv','w')
	f_genotypes = open('genotypes.csv','w')
	f_methylation = open('methylation.csv','w')
	for line in sys.stdin:
		line = line.rstrip()
		header_tokens = []
		f_covariates.write(subjectlist[linenum])
		f_genotypes.write(subjectlist[linenum])
		f_methylation.write(subjectlist[linenum])
		if (linenum==0):
			header_tokens = line.split(delimiter)
			colnum = 0
			for header_token in header_tokens:
				if(colnum>=covariate_start_col and colnum<=covariate_end_col):
					f_covariates.write(delimiter+header_token)
				colnum+=1
			f_covariates.write('\n')
			for snp in snplist:
				f_genotypes.write(','+snp)
			f_genotypes.write('\n')
			for probe in probelist:
				f_methylation.write(','+probe)
			f_methylation.write('\n')
		else:
			print('data line: '+str(linenum))
			tokens = line.split(delimiter)
			colnum = 0
			for token in tokens:
				if(colnum>=covariate_start_col and colnum<=covariate_end_col):
					f_covariates.write(delimiter+token)
					
				elif(colnum==gwas_data_col):
					f_genotypes.write(','+token)
				elif(colnum==meth_data_col):
					f_methylation.write(','+token)
				colnum+=1
			f_covariates.write('\n')
			f_genotypes.write('\n')
			f_methylation.write('\n')
		linenum+=1

if __name__=="__main__":
	main()
