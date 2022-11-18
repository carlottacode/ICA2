#!/usr/bin/python3

import os, sys, subprocess
import pandas as pd
import numpy as np

from functions import enter_protein
from functions import enter_organism
from functions import yes_or_no
from functions import species
from functions import user_input
from functions import choose_from_dict
from functions import seq_prune

import pprint



# =============================================================================
# GET USER INPUT 
# =============================================================================

try:
    
    fseq_headers, prot, orgn = user_input()
    
except:
    print('USER INPUT')


# =============================================================================
# Patmotif
# =============================================================================
try:
    while True:
    
        if yes_or_no('would you like to look for sequence motifs and domains in a \
        sequence you specify in the PROSITE database?')=='y':
            fseq={}
            
            for i in range(0, len(fseq_headers)):
                fseq[i+1]=fseq_headers[i]
            
            pprint.pprint(fseq)
            key, dict_input=choose_from_dict('Enter a number of the sequence you wish to feed into the PROSITE database', fseq)        
    
            temp_file_cmd=f'echo {dict_input}>temp_file'
            pullseq_cmd = f'./pullseq -i file1 -n temp_file > temp_file2'
            remove_line_cmd = f"grep -v '>' temp_file2 > seq_temp"
            patmatmotif_cmd =f"patmatmotifs -sequence seq_temp -outfile seq{key}_patmat.out"
            
            subprocess.call(temp_file_cmd, shell=True)
            subprocess.call(pullseq_cmd, shell=True)
            subprocess.call(remove_line_cmd, shell=True)
            subprocess.call(patmatmotif_cmd, shell=True)
            
            if yes_or_no('Would you like to run another sequence through the PROSITE database?') =='n':
                break
            
        else:
            break
except:
    print('PATMATMOTIF')
   



# =============================================================================
#         #Clustalo
# =============================================================================
try:

    #Our fasta file containing all our sequences is called file1
    #we could use pullseq here to filter out shorter sequences?
    seq_prune()    
    
    if yes_or_no(f'Would you like to align the sequences you retrived from NCBI\
for {prot} in {orgn}?') == 'y':
    
        cmd = "clustalo --threads 12 -i file1 -o clustalo.out"
        subprocess.call(cmd, shell=True)
except:
    print('CLUSTALO')

