import re

sexdict={}
with open("covariates/Demographics_14Dec2025.csv","r") as file:
    header = file.readline()
    #print(f"{header}")
    pattern=r'"(\d+)"'
    patcol=1
    sexcol=9
    for line in file:
        tokens=line.strip().split(",")
        match=re.search(pattern,tokens[patcol])
        patno=match.group(1)
        #print(f"{match.group(1)}")
        if len(tokens[sexcol])==3:
            match=re.search(pattern,tokens[sexcol])
            #print(f"{tokens[sexcol]} {line}")
            sex = int(match.group(1))+1
        else:
            sex = 0
        sexdict[patno] = sex
        #print(f"Adding {sex} to {patno}")


with open("../../../Downloads/PPMI_WGS/pre_imputed/contig_1.psam","r") as file:
    with open("../../../Downloads/PPMI_WGS/pre_imputed/all.psam","w") as out:
        header = file.readline().strip()
        pattern=r'PPMISI(\d+).variant2'
        out.write(f"{header}\n")
        for line in file:
            tokens=line.strip().split("\t")
            match=re.search(pattern,tokens[0])
            if match:
                patno = match.group(1)
                #print(f"patno {patno}")
                if patno in sexdict:
                    out.write(f"PPMISI{patno}.variant2\t{sex}\n")
                else:
                    print(f"ERROR {patno}\t0")
                    out.write(f"PPMISI{patno}.variant2\t0\n")
                    #break
            
