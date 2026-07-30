import sys
import glob
import pandas as pd
import os

direc_in = sys.argv[1]
file_out = sys.argv[2]
print(f"Working on directory {direc_in}")

df_link = pd.read_csv('covariates/ppmi_140_link_list_20210607.csv')
df_link['PATNO'] = df_link['PATNO'].astype('str')

df_cov = pd.read_csv('covariates/Participant_Status_16Jul2025.csv')
df_cov['PATNO'] = df_cov['PATNO'].astype('str')

df_sample_pheno = pd.merge(df_cov,df_link,on='PATNO')[['PATNO','SENTRIXID','POSITION','COHORT_DEFINITION','EVENT_ID']]

#phenos = ["Parkinson's Disease","Healthy Control"]
df_sample_pheno['SENTRIXID'] = df_sample_pheno['SENTRIXID'].astype(str)
df_sample_pheno['SAMPLE_GROUP'] = df_sample_pheno['COHORT_DEFINITION']+'_'+df_sample_pheno['EVENT_ID']
#print(df_sample_pheno)

list_sentrix = []
list_pos = []

sentrix = direc_in.split('/')[-1]
#print(f'direc {direc} sentrix {sentrix}')
#sentrix,pos,_ = fname2.split('_')
#print(f"Filename {sentrix} {pos}")
df_sample_pheno2 = df_sample_pheno.loc[(df_sample_pheno.SENTRIXID==sentrix) ,:]
print(f'df_sample_pheno2\n{df_sample_pheno2}')
if len(df_sample_pheno2) > 0:
    df_sample_pheno2 = df_sample_pheno2.rename(columns={'PATNO':'Sample_Name','SENTRIXID':'Sentrix_ID','POSITION':'Sentrix_Position','SAMPLE_GROUP':'Sample_Group'})
    df_sample_pheno2 = df_sample_pheno2.drop(['COHORT_DEFINITION','EVENT_ID'],axis=1)
    df_sample_pheno2.to_csv(file_out,index=False)

