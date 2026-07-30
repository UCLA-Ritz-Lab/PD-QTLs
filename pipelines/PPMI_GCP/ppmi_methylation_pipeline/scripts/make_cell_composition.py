

mapping_filename='covariates/PPMI_Meth_n524_for_LONI030718.txt'
cell_counts_filename='covariates/EstimateCellCounts_PPMI_EPICn524final030618.csv'
outfilename='covariates/cell_composition.csv'

lookup={}
with open(mapping_filename,'r') as mapping_file:
    header=mapping_file.readline()
    for line in mapping_file:
        tokens = line.rstrip().split('\t')
        lookup[tokens[1]+'_'+tokens[2]] = tokens[0]

with open(cell_counts_filename,'r') as cell_counts_file:
    with open(outfilename,'w') as outfile:
        header=cell_counts_file.readline()
        header=header.replace('Sentrix_position','PATNO')
        outfile.write(header)
        for line in cell_counts_file:
            tokens = line.rstrip().split(',')
            sentrix_pos = tokens[0]
            if sentrix_pos in lookup:
                outfile.write(f"{lookup[sentrix_pos]},{','.join(tokens[1:])}\n")
            
    
            
