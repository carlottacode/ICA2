#!/usr/bin/python3
import os, sys, subprocess
import numpy as np
import pandas as pd

print('''---
Welcome to ICA2
---
''')

#Getting user input so they can specify the proteins and organisms they
#are interested in

try:
	#protein_family=input("What protein family are you interested in?")
	#taxonomic_subset=input("What taxonomic subset do you wish to investigate?")

	protein_family = 'pyruvate dehydrogenase'
	taxonomic_subset='ascomycete fungi'
	esearch_cmd = f'esearch -db protein -query "{protein_family}[PROT]" | grep "Count"'
	line =  subprocess.check_output(esearch_cmd, shell=True).decode("utf-8")

	protein_results = int(line[9:-9])

	if protein_results > 0 :
		print("Luckily for you NCBI has returned "+line[9:-9]+" results for your protein family!")

		#using esearch efilter and efetch to return a fasta file containing sequences
		#in the protein family and the taxonomic subset

		efilter_cmd = f'esearch -db protein -query "{protein_family}[PROT]"|efilter -organism "{taxonomic_subset}"|grep "Count"'
		line =  subprocess.check_output(efilter_cmd, shell=True).decode("utf-8")
		prot_tax_results = int(line[9:-9])
		print(prot_tax_results)
		print(type(prot_tax_results))
		if prot_tax_results > 0:
			print('yes')
			print('This program has found '+str(prot_tax_results)+' results within your specified taxa sub-group')
			print('yes')
			if prot_tax_results > 1000:
				print('''This is a lot of sequences to get from NCBI and so this may take a while
					---
					Make yourself a cup of tea!''')






			cmd =f'esearch -db protein -query "{protein_family}[PROT]"|efilter -organism "{taxonomic_subset}"\
			|efetch -format fasta>file1'
		


			subprocess.check_call(cmd, shell=True)
			try:
				subprocess.call(cmd, shell=True)
				os.system("head file1")
			except:
				print('idk what the issue is')
		else:
			print(protein_family+' returned no results in '+taxonomic_subset+'! Try again!')
	else:
		print("I think you may have spelled your protein family wrong!")
except:
	print("your input wasn't good babe")
finally:
	print('''---
You have reached the end of this program!
---''')
