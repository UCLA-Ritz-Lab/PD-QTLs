#!/usr/bin/python3

import sys
import numpy as np
from sklearn.impute import SimpleImputer
import statsmodels.api as sm
from sklearn.linear_model import LogisticRegression

def main():
	if(len(sys.argv)<4):
		print("Usage: run_epigenome_screen.py [subject disease file] [feature-major file] [outfile]")
		exit(1)
	else:
		diseasefilename=sys.argv[1]
		genomedatafilename=sys.argv[2]
		outdatafilename=sys.argv[3]
		disease_file = open(diseasefilename,'r')
		header = disease_file.readline().rstrip()
		diseaselist = []
		for line in disease_file:
			disease = line.rstrip().split('\t')[1]
			#print('disease',disease)
			diseaselist.append(disease)
		disease_file.close()
		y = np.array(diseaselist).astype('float32')
		#print('disease len',len(diseaselist))
		genome_file = open(genomedatafilename,'r')
		header = genome_file.readline().rstrip()
		imputer = SimpleImputer(strategy='mean')
		model = LogisticRegression()
		fileout=open(outdatafilename,'w')
		for line in genome_file:
			linetokens = line.rstrip().replace('NA','nan').split('\t')
			col1 = linetokens[0]
			#col2 = linetokens[1].replace('NA','nan')
			vals = linetokens[1:]
			#print('val len',len(vals))
			X = np.array(vals).reshape(-1,1)
			X_imputed = imputer.fit_transform(X)
			X2 = sm.add_constant(X_imputed).reshape(-1,2).astype('float32')
			#Xtest = sm.add_constant([[1], [2], [3], [4], [5]])  # Adding a constant for the intercept
			#ytest = [0, 0, 1, 1, 1]
			model = sm.Logit(y,X2)
			result = model.fit()
			#fit = model.fit(X_imputed,y)
			print(col1,'fit',result.summary())
			fileout.write(col1+'\t'+str(result.pvalues[1])+'\n')


		genome_file.close()
		fileout.close()

if __name__=="__main__":
	main()

