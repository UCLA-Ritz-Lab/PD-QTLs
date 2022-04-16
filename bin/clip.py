#!/usr/bin/python3

import sys

def main():
	use_all_fields = False
	user_col_names=set()
	if(len(sys.argv)<3):
		print("Usage: clip.py [column delimiter] [comma separated list of columns]")
		exit(1)
	else:
		delimiter = sys.argv[1]
		if(sys.argv[2] == '*'):
			print("Using all fields",file=sys.stderr)
			print("delimiter:",delimiter," All fields used",file=sys.stderr)
			use_all_fields = True
		else:
			user_col_names_arr=sys.argv[2].split(',')
			print("delimiter:",delimiter," col names:",user_col_names_arr,file=sys.stderr)
			for user_col_name in user_col_names_arr:
				user_col_names.add(user_col_name)
	i=1
	#special_prefix="PDstudy"
	special_prefix=""

	#for user_col_name in user_col_names:
		#print("user col",user_col_name)
	user_col_indices=set()
	for line in sys.stdin:
		line = line.rstrip()
		if(i==1):
			header_tokens = line.split(delimiter)
			#header_tokens = line.split(delimiter)
			col_id = 0
			col_numbers = {}
			printed=False
			for header_token in header_tokens:
				col_numbers[header_token] = col_id
				if(use_all_fields):
					user_col_names.add(header_token)
					
				#print("token:",header_token,col_id,file=sys.stderr)
				#valid_col = False
				#if(header_token.startswith(special_prefix)):
				#	valid_col = True
							#valid_col = True
							#break
				#if(valid_col):
				if (header_token in user_col_names):
					#print("valid col ",header_token,file=sys.stderr)
					if(printed):
						print(delimiter,end='')
					print(header_token,end='')
					user_col_indices.add(col_id)
					#print("col id added ",col_id)
					printed = True
				col_id+=1
			if(printed):
				print()
			#for user_col_index in user_col_indices:
				#print("matched user col index:",user_col_index)
		else:
			line_tokens = line.split(delimiter)
			col_id = 0
			printed=False
			for line_token in line_tokens:
				if(col_id in user_col_indices):
					if(printed):
						print(delimiter,end='')
					if(line_token == '' or line_token == 'NA'):
						line_token = '\\N'
					elif(line_token == '.'):
						line_token = '0.0'
					print(line_token,end='')
					#print("In set: ",col_id,line_token,end='')
					printed = True
				col_id+=1
			if(printed):
				print()
			#print(line)
		i+=1

if __name__=="__main__":
	main()
