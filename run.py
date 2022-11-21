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
from functions import define_subset
import pprint
import textwrap
from functions import print_info
from functions import print_graph
from functions import print_process
from functions import print_checkpoint 
from functions import clustalo
from functions import hmoment


subprocess.call('mkdir prog_run', shell=True)
output_path = './prog_run/'


if subprocess.check_output('cat ~/.bash_profile | grep edirect', shell=True).decode("utf-8") == '':
    print_info('Hi there, it looks like you don\'t have Entrez Direct installed which is required to run this program...')
    if yes_or_no('Would you like to install Entrez Direct?') == 'y':
        subprocess.call('sh -c "$(wget -q ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)"', shell=True)
    


# =============================================================================
# GET USER INPUT 
# =============================================================================
try:
    
    headers, prot, orgn = user_input()
    
    fseq={}
    for i in range(0, len(headers)):
        fseq[i+1]=headers[i]
    # =============================================================================
    #         #Clustalo
    # =============================================================================
    try:
        clustalo(prot, orgn)
        
        # =============================================================================
        #         #Infoalign
        # =============================================================================
        try:
            cmd = f"infoalign -noname -sequence {output_path}clustalo.phy -outfile {output_path}infoalign_out.txt"
            subprocess.call(cmd, shell=True)
            
            
            with open(f'{output_path}infoalign_out.txt') as myfile:
                    infoalign_df=pd.read_csv(myfile, sep='\t')
                    infoalign_df.sort_values(['% Change'], inplace=True)
                    infoalign_df.to_csv(f'{output_path}infoalign_out.csv')
            
            if yes_or_no('Would you like to see some information about the alignment?') == 'y':
                    print(infoalign_df.head())
                    my_min = infoalign_df['% Change'].loc[infoalign_df['% Change'].idxmin()]
                    
                    if my_min > 50:
                        print_info(f'It looks like your sequences are all more than >{my_min}< % different from eachother. Don\'t fret this may just indicate that these proteins are not very conserved. Let\'s have a look at a plot. Let\'s visualise this in a plot!)')
                    else:
                        print_info(f'It looks like your sequences are all at least {100-my_min} similar to eachother! Let\'s visualise this in a plot!')
                        
            # =============================================================================
            #         #Plotcon
            # =============================================================================
            try:
                if yes_or_no('Would you like to continue?') =='y':
                    plotcon_create_cmd = f"plotcon -winsize 6 -sequence {output_path}clustalo.phy -auto -graph png -stdout"
                    plotcon_view_cmd = f"plotcon -winsize 6 -sequence {output_path}clustalo.phy -auto -graph x11"
                    subprocess.call(plotcon_create_cmd, shell=True)
                    subprocess.call(plotcon_view_cmd, shell=True)
                    
                    print_info('If there were regions in the protein sequences which appeared to be highly conserved there may be specific domains or motifs present in this group of proteins... Let\'s explore this further.')
                    
                    # =============================================================================
                    # Patmotif
                    # =============================================================================
                    try:
                        if yes_or_no('Would you like to use the PROSITE database to search for motifs and domains in your sequences?')=='y':
                            if yes_or_no('Would you like to look at motifs in the consensus sequence?')=='y':
                                
                                cons_cmd = f'cons {output_path}clustalo.phy {output_path}cons_seq'
                                subprocess.call(cons_cmd, shell=True)
                                pat_cmd = f"patmatmotifs {output_path}cons_seq {output_path}{(prot[:5]+'_'+orgn[:5]).replace(' ', '')}_con.motif"
                                subprocess.call(pat_cmd, shell=True)
                            
                            while True:
                                if yes_or_no('Would you like to run a different sequence through the PROSITE database?') =='y':
                                    
                                    
                                    pprint.pprint(fseq)
                                    key, dict_input=choose_from_dict('Enter a number of the sequence you wish to feed into the PROSITE database...', fseq)        
                            
                                    temp_file_cmd=f'echo {dict_input}>{output_path}temp_file'
                                    pullseq_cmd = f'./pullseq -i {output_path}file1 -n {output_path}temp_file > {output_path}temp_file2'
                                    remove_line_cmd = f"grep -v '>' {output_path}temp_file2 > {output_path}seq_temp"
                                    patmatmotif_cmd =f"patmatmotifs -sequence {output_path}seq_temp -outfile {output_path}seq_{key}.motif"
                                    
                                    subprocess.call(temp_file_cmd, shell=True)
                                    subprocess.call(pullseq_cmd, shell=True)
                                    subprocess.call(remove_line_cmd, shell=True)
                                    subprocess.call(patmatmotif_cmd, shell=True)
                                else:
                                    break
                        
                        # =============================================================================
                        # Phylogenetic tree
                        # =============================================================================
                        try:
                            if yes_or_no('Would you like to see the tree of proetin sequences you parsed into the multiple sequence aligner?') =='y':
                                accn = []
                                full_header=[]
                                
                                for i in headers:
                                    i = i.rstrip()
                                    
                                    headers_F = i.replace('[', '').replace(']', '').replace(' ', '_').replace(':', '')
                                    full_header.append(headers_F)
                                    
                                    i = i.split(' ')
                                    accn.append(i[0])
                                
                                accession = dict(zip(accn, full_header))
                            
                                with open(f'{output_path}clustalo.tree') as myfile:
                                    tree = myfile.read()
                                    tree = tree.replace('\n','')
                                    for i in accn:
                                        tree = tree.replace(i, accession[i])
                
                                with open(f'{output_path}outtree', 'w') as myfile:
                                    myfile.write(tree)
                                    
                                cmd_create_tree = f"figtree -graphic PDF {output_path}outtree {output_path}outtree.pdf"
                                cmd_view_tree = f'ghostscript {output_path}outtree.pdf'
                                subprocess.call(cmd_create_tree, shell=True)
                                subprocess.call(cmd_view_tree, shell=True)
                                
                            # =============================================================================
                            # HMOMENT      
                            # =============================================================================
                            try:
                                hmoment(fseq)
                            except:
                                print('Something went wrong trying to use hmoment... Please try again')
                        except:
                            print('Phylogenetic tree didn\'t work')    
                        
                    except:
                        print('Something went wrong trying to use patmatmotifs... Please try again')
            except:
                print('Something went wrong trying to use plotcon... Please try again')    
            
        except:
            print('Something went wrong trying to use infoalign... Please try again')     
    
                
    except:
        print('Seomthing went wrong with using clustalo... Please try again')
    
except:
    print('Something went wrong with the user input... Please try again')






   


        
