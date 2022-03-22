## Importing rawdata into SQL tables

First obtain files from https://app.box.com/folder/158643387254

### STEP 1: Loading methylation Rdata objects in

Download the files under the Data folder.  Unzip files so that the contents are under <repo_root>/rawdata. In <repo_root>/bin, edit the files rdata2tsv.sh and tsv2sql.sh so that the project_dir variable points to the <repo_root>.

Under <repo_root>/rawdata/Methylation, for each of the RData files, you can run:

 ../../bin/rdata2tsv.sh

which will give you the arguments necessary for converting the RData objects into text files that the SQL database can parse on import.  Note that the first argument for this shell script is simply the RData file name without the .RObject extension.  The second argument is the key type for the rows of the RObject and the third argument is the key type for the columns.  So if rows corresponded to subjects and columns as probes, you can specify 'subject' for the second argument and 'probe' for the third argument.  Once you have converted each of the RObjects, you will see the same number of new files as original RObject files in the same directory.  You can then proceed towards importing these into SQL tables.


Under <repo_root>/rawdata/Methylation, run 
 
 ../../bin/tsv2sql.sh

which will provide necessary arguments for this shell script.  The first argument is simply the name of the new tsv file that was generated in the previous step.  The second column will be the index that represents the rows in the matrix in the tsv file.  So the second argument could be subject_id or probe_id for example.

A convenience script that loops over all the tsv files in the current directory is provided in <repo_root>/bin/ called load_all_tsv.sh.  This script can just be simply run in the same folder as the tsv files are located in, without any arguments passed in.  

### STEP 2: Loading methylation annotation in

Under <repo_root>/rawdata/Methylation, run the Python script
 
 ../../bin/clip.py

which will provide arguments to selecting the necessary columns from a field terminated text file (in this case the Illumina annotation file) into a format that is friendly for MySQL import (e.g. Null fields will be printed as \N). In this case you can run:

  ../../bin/clip.py , '*' < HumanMethylation450_15017482_v.1.2.csv > illumina_annotation_import.csv

which would save the new file as illumina_annotation_import.csv in the current folder.

There is a convenience Python script that generates a SQL import script, which in turn needs to be slightly edited so that certain db table columns can have the correct data type. This can be run as:

 ../../bin/make_create_table.py

to display help on arguments.  In our case we can run:

 ../../bin/make_create_table.py illumina_annotation ',' < illumina_annotation_import.csv  > ../../bin/create_illumnina_annotation.sql

which will save a SQL script in the bin folder. Open this script with vi to edit the file accordingly :

 vi ../../bin/create_illumnina_annotation.sql 

Certain things to edit: Columns such as numerics can be converted from text for float for example.  Also, make sure the load data infile statement points to the illumina_annotation_import.csv file. We can create indices too. Make sure CHR is changed to chrom, and its datatype is varchar(10).  MAPINFO should be a an int type.  The string ", unique key index_map(chrom,MAPINFO)" should be added right after the string "primary key ( IlmnID ),".  Now run the script as:

 sql_pd_qtl < ../../bin/create_illumnina_annotation.sql




