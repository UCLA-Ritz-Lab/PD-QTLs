import pandas as pd
import glob
import sys
import os

idat_out=snakemake.input.idat_output

#matrices=['beta','M','pd','intensity','detP','beadcount','Meth','UnMeth']
matrices=['beta']
samplesheets = ['pd']

def merge_files(pd_file,rows_file,vals_filelist,outfilename):
    fileout = open(outfilename,'w')    
    with open(pd_file,'r') as file:
        file.readline()
        colnames = [line.strip().split('\t')[0] for line in file.readlines()]
        fileout.write('\t'.join(colnames))
    fileout.write('\n')    
    
    with open(rows_file,'r') as file:
        rownames = [line.strip() for line in file.readlines()]    
    filehandles = []
    for filename in vals_filelist:
        filehandles.append(open(filename,'r'))
    for index,rowname in enumerate(rownames):
        fileout.write(rowname)
        totalcols=0
        for index2,filehandle in enumerate(filehandles):
            fileout.write('\t')
            linein = filehandle.readline().strip()
            #print(f"index {index2} linein {linein}")
            cols = len(linein.split('\t'))
            totalcols+=cols
            fileout.write(linein)
        fileout.write('\n')
        #print(f"total cols {totalcols}")
        if index % 5000 == 0:
            print(f"Row {index} completed")        
    for filehandle in filehandles:
        filehandle.close()
    fileout.close()

def merge_samplesheets(cols_file,vals_filelist,outfilename):
    with open(cols_file,'r') as file:
        colnames = [line.strip() for line in file.readlines()]
    fileout = open(outfilename,'w')    
    for index,colname in enumerate(colnames):
        if index>0:
            fileout.write('\t')
        fileout.write(colname)
    fileout.write('\n')
    for filename in vals_filelist:
        with open(filename,'r') as file:
            for line in file:
                linetokens = line.strip().split('\t')
                samplegroupvisit = linetokens[3]
                samplegroup,visit = samplegroupvisit.split('_')                                
                samplegroup = samplegroup.replace("'","").replace(" ",".")
                fileout.write(linetokens[0]+'_'+visit+'\t'+linetokens[1]+'\t'+linetokens[2]+'\t'+samplegroup+'.'+visit)
                fileout.write('\n')
    fileout.close()
    
        
    
sentrixlist = glob.glob(idat_out+"/*")
def process_datasets(datasets,orientation):
    for matrix in datasets:    
        cols_filelist=[]
        vals_filelist=[]
        for index,sentrix in enumerate(sentrixlist):
            if index==0:                
                    rows_file=sentrix+'/'+matrix+'_row.txt'
            #print(f"sentrix {index} {sentrix}")            
            cols_file=sentrix+'/'+matrix+'_cols.txt'
            cols_filelist.append(cols_file)
            vals_file=sentrix+'/'+matrix+'_vals.txt'
            vals_filelist.append(vals_file)
        outfilename=matrix+'.txt'
        print(f"Output to {outfilename}")
        #print(f"{cols_file} {rows_file} {outfilename} {vals_filelist} ")
        if orientation == 'vertical':
            merge_samplesheets(cols_file,vals_filelist,outfilename)
        elif orientation == 'horizontal':
            merge_files('pd.txt',rows_file,vals_filelist, outfilename)
            
            
process_datasets(samplesheets,'vertical')
process_datasets(matrices,'horizontal')
