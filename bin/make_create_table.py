#!/usr/bin/python3

import sys

dbname = 'pd_qtl'
data_dir = '/home/garyc/analysis/PD-QTLs/rawdata'

def main():
	if(len(sys.argv)<3):
		print("Summary: this convenience script generate a template sql script for importing csv files that the user can further edit.  Be sure to edit set the dbname and data_dir variable at the top of this script.\nUsage: make_create_table.py [table_name] [delimiter]")
		exit(1)
	else:
		table_name = sys.argv[1]
		delimiter = sys.argv[2]
	#for user_col_name in user_col_names:
		#print("user col",user_col_name)
	user_col_indices=set()
	print("use "+dbname+";");
	print("drop table if exists ",table_name,";")
	print("create table ",table_name," (",end='')
	for line in sys.stdin:
		line = line.rstrip()
		break;
	table_cols = line.split(delimiter)
	i=1
	for table_col in table_cols:
		if(i==1):
			pk = table_col
			print(table_col.replace(".","_")," varchar(255), ",end='')
		else:
			print(table_col.replace(".","_")," text, ",end='')
		i+=1
	print(" primary key (",pk,"));")
	print("load data infile '"+data_dir+"' into table ",table_name," fields terminated by ',' ignore 1 lines;");

if __name__=="__main__":
	main()

