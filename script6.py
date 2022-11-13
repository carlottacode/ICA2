#!/usr/bin/python3

import os, sys, subprocess
import numpy as np
import pandas as pd

# =============================================================================
# YES OR NO FUNCTION
# =============================================================================
def yes_or_no(question):
    counter = ''
    while counter == '':
        
        user = input(question+'(y/n)')
        
        if user.lower() in ['yes', 'y']:
            counter += 'y'
            return counter
        elif user.lower() in ['no', 'n']:
            counter+='n'
            return counter
        else:
            print('I\'m sorry I don\'t understand...')
            print('Please enter y/n in the accepted format')
            
# =============================================================================
# PROTEIN FUNCTION
# =============================================================================
def enter_protein():
    while True:
        protein_family=input("What protein family are you interested in?")
        #protein_family = 'pyruvate dehydrogenase'
        
        #Checking that the protein family returns results from NCBI if not then 
        #there may be spelling mistakes 
        esearch_cmd = f'esearch -db protein -query "{protein_family}[PROT]" | grep "Count"'
        line =  subprocess.check_output(esearch_cmd, shell=True).decode("utf-8")
        protein_results = int(line[9:-9])
        
        if protein_results == 0:
            print('''---
                  Hi there! Your protein didn\'t return any results from NCBI...
                  This means you might have made a typo...
                  Please try again!
                  ---''')
        elif 1<= protein_results <= 2:
            print(f'---\nHi ! NCBI returned {protein_results} which isn\'t enough for the subsequent analysis.\
                  I\'m sorry about this...\n Please try again!')
                
        elif protein_results > 2 :
                print(f"Luckily for you NCBI has returned {line[9:-9]} sequence results for your protein family {protein_family}!")
    
                #using esearch efilter and efetch to return a fasta file containing sequences
                #in the protein family and the taxonomic subset
                break
        else:
            print("Something isn't working please try again...")
    
    return protein_family


# =============================================================================
# ORGANISM FUNCTION
# =============================================================================
def enter_organism(protein):
    while True:
        taxonomic_subset=input("What taxonomic subset do you wish to investigate?")
        #taxonomic_subset='ascomycete fungi'
        
        efilter_cmd = f'esearch -db protein -query "{protein}[PROT]"|efilter -organism "{taxonomic_subset}"|grep "Count"'
        
        line =  subprocess.check_output(efilter_cmd, shell=True).decode("utf-8")
        prot_tax_results = int(line[9:-9])
        
        if prot_tax_results == 0:
            print(f'''---
                  Hello,
                  Your combination of >{taxonomic_subset}< and >{protein}< didn\'t return any results from NCBI...
                  This means you might have made a typo...
                  Please try entering the taxonomic subgroup again!
                  ---''')
        elif 1 <= prot_tax_results <= 2:
            print(f'---\nHi ! NCBI returned {line[9:-9]} which isn\'t enough for the subsequent analysis.\
                  I\'m sorry about this...\n Please try again!')
            
        elif 2 < prot_tax_results < 1000:
                print(f'This program has found {line[9:-9]} results within your specified taxa sub-group\
                      We can now continue with the analysis...')
                break
                
        elif prot_tax_results >= 1000:
            print('This is a lot of sequences to get from NCBI and so this may take a while...')
            if yes_or_no('Would you still like to continue?') == 'y':
                print('Make yourself a cup of tea!')
                break
            elif yes_or_no('Would you like to enter a different protein family?')=='y':
                protein = enter_protein()
        
        else:
            print("Something is going wrong, please try again!")
                                    
    return protein, taxonomic_subset

#Getting user input so they can specify the proteins and organisms they
#are interested in

try:
# =============================================================================
#   #Error trapping protein family user input...    
# =============================================================================
    try:    
        prot=enter_protein()
    except:
        print('PROTEIN FUNCTION')
# =============================================================================
#   #Error trapping taxonomic subset user input...
# =============================================================================
    try:    
        prot2, orgn=enter_organism(prot)    
        print(prot2) 
        print(orgn)
    except:
        print('ORGANISM FUNCTION')
        
except:
    print('WHOLE THING')
