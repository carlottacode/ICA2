#!/usr/bin/python3
import os, sys, subprocess
import numpy as np
import pandas as pd

print('''Welcome to ICA2
---
''')

#Getting user input so they can specify the proteins and organisms they
#are interested in

try:
	protein_family=input("What protein family are you interested in?")
	taxonomic_subset=input("What taxonomic subset do you wish to investigate?")

	#protein_family = 'pyruvate dehydrogenase'
	#taxonomic_subset='ascomycete fungi'
	esearch_cmd = f'esearch -db protein -query "{protein_family}[PROT]" | grep "Count"'
	line =  subprocess.check_output(esearch_cmd, shell=True).decode("utf-8")


	if int(line[9:-9]) > 0 :


		#using esearch efilter and efetch to return a fasta file containing sequences
		#in the protein family and the taxonomic subset
		cmd =f'esearch -db protein -query "{protein_family}[PROT]"|efilter -organism "{taxonomic_subset}"\
		|efetch -format fasta>file1'

		subprocess.check_call(cmd, shell=True)
		try:
			subprocess.call(cmd, shell=True)
			os.system("head file1")
		except:
			print('idk what the issue is')
except:
	print("your input wasn't good babe")
finally:
	print("That's it!")
