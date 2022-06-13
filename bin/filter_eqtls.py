#!/usr/bin/python3

import sys

def get_info_set(list_file):
	f_handle = open(list_file,'r')
	snp_set = set()
	for line in f_handle:
		line = line.strip()
		rs_allele = line.split('_')
		rs = rs_allele[0]	
		snp_set.add(rs)
	f_handle.close()
	return snp_set

def print_info_list(info_list):
	for info in info_list:
		print(info)

def main():
	# some constants
	delimiter = '\t'
	if(len(sys.argv)<3):
		print("Usage: "+sys.argv[0]+" [whitelist|blacklist] [noheader snpfile]")
		exit(1)
	list_mode = sys.argv[1]
	list_file = sys.argv[2]
	snp_set = get_info_set(list_file)
	linenum = 0
	for line in sys.stdin:
		line = line.rstrip()
		header_tokens = []
		if (linenum==0):
			print(line)
		elif (linenum>0):
			tokens = line.split(delimiter)
			colnum = 0
			for token in tokens:
				if(colnum==1):
					if((list_mode=='whitelist' and token in snp_set) or (list_mode=='blacklist' and token not in set_set)) :
						print(line)
						
				colnum+=1
		linenum+=1

if __name__=="__main__":
	main()
