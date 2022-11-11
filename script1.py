#!/usr/bin/python3
import os, sys, subprocess
import numpy as np
import pandas as pd

print('''Welcome to ICA2
---
''')

#Getting user input so they can specify the proteins and organisms they
#are interested in

#protein_family=raw_input("What protein family are you interested in?")
#taxonomic_subset=raw_input("What taxonomic subset do you wish to investigate?")

protein_family = 'pyruvate dehydrogenase'
taxonomic_subset='ascomycete fungi'


#using esearch efilter and efetch to return a fasta file containing sequences
#in the protein family and the taxonomic subset
cmd =f'esearch -db protein -query "{protein_family}[PROT]"|efilter -organism "{taxonomic_subset}"\
|efetch -format fasta>file1'

subprocess.call(cmd, shell=True)

os.system("head file1")
