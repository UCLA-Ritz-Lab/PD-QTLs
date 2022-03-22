Importing rawdata into SQL tables
================

## STEP 1: Obtain files from https://app.box.com/folder/158643387254

Download the files under the Data folder.  Unzip files so that the contents are under <repo_root>/rawdata. In <repo_root>/bin, edit the files rdata2tsv.sh and tsv2sql.sh so that the project_dir variable points to the <repo_root>.

Under <repo_root>/rawdata/Methylation, for each of the RData files, you can run:

 ../../bin/rdata2tsv.sh

which will give you the arguments necessary for converting the RData objects into text files that the SQL database can parse on import.  Note that the first argument for this shell script is simply the RData file name without the .RObject extension.  The second argument is the key type for the rows of the RObject and the third argument is the key type for the columns.  So if rows corresponded to subjects and columns as probes, you can specify 'subject' for the second argument and 'probe' for the third argument.  Once you have converted each of the RObjects, you will see the same number of new files as original RObject files in the same directory.  You can then proceed towards importing these into SQL tables.


Under <repo_root>/rawdata/Methylation, run 
 
 ../../bin/tsv2sql.sh

which will provide necessary arguments for this shell script.  The first argument is simply the name of the new tsv file that was generated in the previous step.  The second column will be the index that represents the rows in the matrix in the tsv file.  So the second argument could be subject_id or probe_id for example.

A convenience script that loops over all the tsv files in the current directory is provided in <repo_root>/bin/ called load_all_tsv.sh.  This script can just be simply run in the same folder as the tsv files are located in, without any arguments passed in.  
