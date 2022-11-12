#!/usr/bin/python3

import pandas as pd
import numpy as np
import os, sys, subprocess
from script1 import protein_family, taxonomic_subset

print('''---
We're going to start some analysis of your sequences with EMBOSS programs...
---''')

#Clustalo
#Our fasta file containing all our sequences is called file1
#we could use pullseq here to filter out shorter sequences?


cmd = "clustalo --threads 12 -i file1 -o clustalo.out"
subprocess.call(cmd, shell=True)

#Infoalign

cmd = "infoalign -noweight -noname -sequence clustalo.out -outfile infoalign_out.txt"
subprocess.call(cmd, shell=True)

with open('./infoalign_out.txt') as myfile:
	infoalign_df=pd.read_csv(myfile, sep='\t')

print(infoalign_df.head())

#Showalign
output_file = protein_family[:5]+taxonomic_subset[:5]
cmd =f"showalign -order S -show S -sequence clustalo.out -outfile {output_file}.showalign| cat -"
subprocess.call(cmd, shell=True)

#Plotcon
cmd = "plotcon -winsize 4 -sequence clustalo.out -graph svg"

#Patmotif
cmd = "patmatmotif "
