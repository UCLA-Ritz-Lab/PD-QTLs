#!/usr/bin/python3

import sys

def main():
	delimiter = '\t'
	comma = ','
	header = False
	f_snpinfo = open('snpinfo.txt','w')
	f_genotypes = open('genotypes.txt','w')
	for line in sys.stdin:
		line = line.rstrip()
		tokens = line.split(delimiter)
		if (len(tokens)>1):
 			#header found
			if(tokens[0] == '#CHROM'):
				header = True
				f_subjects = open('subjects.txt','w')
				i = 0
				for token in tokens:
					if(i>9):
						f_subjects.write(token+'\n')
					i+=1
				f_subjects.close()
			else:
				if(header):
					#this should be the body
					f_snpinfo.write(tokens[2]+'\t'+tokens[0]+'\t'+tokens[1]+'\t'+tokens[3]+'\t'+tokens[4]+'\n')
					f_genotypes.write(tokens[2]+delimiter+comma.join(tokens[9:len(tokens)])+'\n')
					
	f_snpinfo.close()
	f_genotypes.close()

if __name__=="__main__":
	main()
