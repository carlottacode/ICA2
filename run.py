#!/usr/bin/python3

import os, sys, subprocess  
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from functions import enter_protein
from functions import enter_organism
from functions import yes_or_no
from functions import species
from functions import user_input
from functions import choose_from_dict
from functions import prune_seq
from functions import check_num

import pprint

subprocess.call('mkdir prog_run', shell=True)
output_path = './prog_run/'



# =============================================================================
# GET USER INPUT 
# =============================================================================

try:
    
    fseq_headers, prot, orgn = user_input()
    
except:
    print('Something went wrong with the user input... Please try again :)')



# =============================================================================
#         #Clustalo
# =============================================================================
try:

    #Our fasta file containing all our sequences is called file1
    #we could use pullseq here to filter out shorter sequences?  
    if yes_or_no(f'Would you like to align the sequences you retrived from NCBI\
for {prot} in {orgn}?') == 'y':
    
        
        with open(f'{output_path}file_wo_partial') as myfile:
            
            seq_lines = [line for line in myfile]
        
        seq =[]  
        for line in seq_lines:
            if '>' in line:
                seq.append(line.replace('\n', '>'))
            else:
                seq.append(line)
                
            
        
        seq=''.join(seq)
        
        
        seq = seq.replace('X', '')
        seq = seq.replace('\n', '')
        
        seq = seq.split('>')
        
        seq = seq[1:]

        seq_df = pd.DataFrame(seq[1::2], index=seq[::2], columns=['Sequence'])
        seq_df['Sequence Length'] = seq_df['Sequence'].str.len()
        
        my_max = seq_df['Sequence Length'].loc[seq_df['Sequence Length'].idxmax()]      # Maximum in column
        my_min = seq_df['Sequence Length'].loc[seq_df['Sequence Length'].idxmin()]      # Maximum in column
        
        seq_df.sort_values(['Sequence Length'], ascending=False, inplace=True)
        
        print(seq_df)
        
        plt.hist(seq_df['Sequence Length'])
        plt.title('Sequence Lengths for your Search Query')
        plt.savefig(output_path+'seq_histogram.pdf', format='pdf')
        
        cmd_view_hist = f'gs {output_path}seq_histogram.pdf'
        subprocess.call(cmd_view_hist, shell=True)
        
        while True:
            num_seq = len(seq_df.index.tolist())
            if yes_or_no('Would you like to specify a subset?') == 'n':
                print(f'\nContinuing with all {num_seq} sequences...\n'.center(150, '-'))
                break
            else: 
            
                if yes_or_no('Would you like to remove some of the shorter sequences?') == 'y':
                    num = check_num(my_max, my_min, 'Please enter a minimum sequence length.')
                
                    
                    
                    subset=seq_df[seq_df['Sequence Length']>num]
                    num_seq = len(subset.index.tolist())
                    
                    seq = subset.index.tolist()
                    if num_seq > 2:
                        print('\nContinuing with this subset!\n'.center(150,'~'))
                        break
                    else:
                        print(f'\nI\'m really sorry but your subset only contains\
{num_seq} sequences and you need at least three to to a multiple sequence alignemnt\n'.center(150, '-'))
                
                    print(f'With this constraint you want to align {num_seq} sequences')
                if yes_or_no('Would you like to remove some of the longer sequences?') == 'y':
                    num = check_num(my_max, my_min, 'Please enter a maximum sequence length.')
                
                    #seq_df.loc['KAG0683330.1 helicase [[Candida]']
                    subset=seq_df[seq_df['Sequence Length']<num]
                    num_seq = len(subset.index.tolist())
                    seq = subset.index.tolist()
                    print(f'With this constraint you want to align {num_seq} sequences')
                    if num_seq > 2:
                        print('\nContinuing with this subset!\n'.center(150,'~'))
                        break
                    else:
                        print(f'\nI\'m really sorry but your subset only contains\
{num_seq} sequences and you need at least three to to a multiple sequence alignemnt\n'.center(150, '-'))
                    
                
                    print(f'With this constraint you want to align {num_seq} sequences')
                if yes_or_no('Would you like to define your own range?') == 'y':
                    min_num = check_num(my_max, my_min, 'Please enter your minimum sequence length.')
                    max_num = check_num(my_max, my_min, 'Please enter your maximum sequence length.')
                    subset=seq_df[seq_df['Sequence Length'].between(min_num, max_num)]
                    num_seq = len(subset.index.tolist())
                    seq = subset.index.tolist()
                    print(f'With this constraint you want to align {num_seq} sequences')
                    if num_seq > 2:
                        print('\nContinuing with this subset!\n'.center(150,'~'))
                        break
                    else:
                        print(f'\nI\'m really sorry but your subset only contains\
{num_seq} sequences and you need at least three to to a multiple sequence alignemnt\n'.center(150, '-'))
                    
        seq_file      = []
        for i in seq:
            i = i+'\n'
            seq_file.append(i)               
            
        with open(f'{output_path}msa_in', 'w') as myfile:
            myfile.writelines(seq_file)
        
        pullseq_cmd = f'./pullseq -i {output_path}file1 -n {output_path}msa_in>{output_path}clustalo_in'  
        subprocess.call(pullseq_cmd, shell = True)
        cmd = f"clustalo --force --full --threads 16 --outfmt=phy \
            --guidetree-out={output_path}clustalo.tree -i {output_path}clustalo_in -o {output_path}clustalo.phy"

        subprocess.call(cmd, shell=True)
except:
    print('CLUSTALO')
    
# =============================================================================
# Phylogenetic tree
# =============================================================================
        

try:
    
    if yes_or_no('Would you like to see the tree of proetin sequences you parsed into\
the multiple sequence aligner?') =='y':
    
        with open(f"{output_path}headers_wo_partial.txt") as myfile:
            lines = [line for line in myfile]
            
        
        accn = []
        header=[]
        
        for i in lines:
            i = i.rstrip()
            
            headers_F = i.replace('[', '')
            headers_F = headers_F.replace(']', '').replace(' ', '_').replace(':', '')
            header.append(headers_F)
            
            i = i.split(' ')
            accn.append(i[0])
        
        
        accession = dict(zip(accn, header))
        pprint.pprint(accession)
    
        with open(f'{output_path}clustalo.tree') as myfile:
            tree = myfile.read()
        
        tree = tree.replace('\n','')
        print(tree)
        for i in accn:
            tree = tree.replace(i, accession[i])
            
        with open(f'{output_path}outtree', 'w') as myfile:
            myfile.write(tree)
            
        cmd_create_tree = f"figtree -graphic PDF {output_path}outtree {output_path}outtree.pdf"
        cmd_view_tree = f'ghostscript {output_path}outtree.pdf'
        subprocess.call(cmd_create_tree, shell=True)
        subprocess.call(cmd_view_tree, shell=True)
except:
    print('Phylogenetic tree didn\'t work')    



# =============================================================================
#         #Infoalign
# =============================================================================
        
try:
    if yes_or_no('Would you like to see some information about the alignment?') == 'y':
            
            cmd = f"infoalign -noname -sequence {output_path}clustalo.phy -outfile {output_path}infoalign_out.txt"
            subprocess.call(cmd, shell=True)
            
            with open(f'{output_path}infoalign_out.txt') as myfile:
                    infoalign_df=pd.read_csv(myfile, sep='\t')
            
            print(infoalign_df.head())
except:
    print('INFOALIGN')        
    
# =============================================================================
#         #Showalign
# =============================================================================
try:
    
        if yes_or_no('Would you like to visualise this alignment?')=='y':
            output_file = prot[:5]+'_'+orgn[:5]
            output_file.replace(' ', '')
            cmd =f"showalign -order S -show S -sequence {output_path}clustalo.phy -outfile {output_path}{output_file}.showalign| cat -"
            subprocess.call(cmd, shell=True)
except:
    print('Showalign didn\'t work')
        
# =============================================================================
#         #Plotcon
# =============================================================================
try:
        if yes_or_no("Would you like to see a plot of the conservation between your protein sequences?")=='y':
            
            cmd = f"plotcon -winsize 6 -sequence {output_path}clustalo.phy -auto -graph x11"
            subprocess.call(cmd, shell=True)
except:
    print('Plotcon didn\'t work')    
               
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
    
            temp_file_cmd=f'echo {dict_input}>{output_path}temp_file'
            pullseq_cmd = f'./pullseq -i file1 -n {output_path}temp_file > {output_path}temp_file2'
            remove_line_cmd = f"grep -v '>' {output_path}temp_file2 > {output_path}seq_temp"
            patmatmotif_cmd =f"patmatmotifs -sequence {output_path}seq_temp -outfile {output_path}seq{key}_patmat.out"
            
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

