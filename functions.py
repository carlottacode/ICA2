#!/usr/bin/python3
import os, sys, subprocess
import numpy as np
import pandas as pd

print('\n---\nWelcome to ICA2\n---\n')


# =============================================================================
# This function allows you to ask a yes or no question which either returns a 
# counter 'y' or 'n' which you can make use of in other functions when deciding 
# what to do next i.e. whether you want to exit a while True: loop
# =============================================================================
def yes_or_no(question):
    counter = ''
    while counter == '':
        
        user = input('\n\n!!!\n'+question+'(y/n)\n!!!\n')
        
        if user.lower() in ['yes', 'y']:
            counter += 'y'
            return counter
        elif user.lower() in ['no', 'n']:
            counter+='n'
            return counter
        else:
            print('\n---\nI\'m sorry I don\'t understand... Please enter y/n \
in the accepted format\n---\n')
            

# =============================================================================
# This function asks the user to enter a protein family of interest. It will 
# not allow the user to continue unless the protein family returns sequences
# from NCBI. Therefore hopefully catching all user input errors.
# =============================================================================
def enter_protein():
        while True:
            
                protein_family=input("\n###\nWhat protein family are you interested in?\n###\n").lower()
                #protein_family = 'pyruvate dehydrogenase'
                
                if protein_family[-1] == 's':
                  protein_family = protein_family[:-1]
                  
                #Checking that the protein family returns results from NCBI if not then 
                #there may be spelling mistakes 
                esearch_cmd = f'esearch -db IPG -query \
                    "{protein_family}[PROT]" | grep "Count"'
                
                line =  subprocess.check_output(esearch_cmd, shell=True).decode("utf-8")
                protein_results = int(line[9:-9])
                
                if protein_results == 0:
                    print('\n')
                    print(f"\nHi there! Your protein (>{protein_family}<) \
didn\'t return any results from NCBI...\nThis means you might have made a typo... \n\
Please try again!\n".center(350, '-'))
                              
                elif 1<= protein_results <= 2:
                    print('\n')
                    print(f'\nHi ! NCBI returned >{protein_results}< which \
isn\'t enough for the subsequent analysis. I\'m sorry about this... \n\
Please try again!\n'.center(350, '-'))
                        
                elif protein_results > 2 :
                    print('\n')
                    print(f"\nLuckily for you NCBI has returned >{line[9:-9]}< \
sequence results for your protein family (>{protein_family}<)!\n".center(300, '-'))
    
                    break
                else:
                    print("Something isn't working please try again...")
        
        return protein_family


# =============================================================================
# This functions asks the user to enter a taxonomic subset. It doesn't allow 
# the user to continue with the analysis until a combination of protein family 
# and taxonomic subset returns more than 2 sequences from NCBI. There are points
# at which if the number of sequences returned is not satisfactory that you are
# able to enter a different protein family if you wish, starting again...
# =============================================================================
def enter_organism(protein):
        while True:
                taxonomic_subset=input("\n###\nWhat taxonomic subset do you wish to investigate?\n###\n").lower()
                #taxonomic_subset='ascomycete fungi'
                
                efilter_cmd = f'esearch -db IPG -query "{protein}[PROT]\
                    AND {taxonomic_subset}[ORGN]"|grep "Count"'
                
                line =  subprocess.check_output(efilter_cmd, shell=True).decode("utf-8")
                prot_tax_results = int(line[9:-9])
                
                if prot_tax_results == 0:
                    print('\n')
                    print(f"\nHello, your combination of >{taxonomic_subset}< \
and >{protein}< didn\'t return any results from NCBI...\nThis means you might \
have made a typo...\nPlease try entering the taxonomic subgroup again!\n".center(350, '-'))
                    
                    if yes_or_no('Would you like to enter a different protein family?') == 'y':
                        protein = enter_protein()
                              
                elif 1 <= prot_tax_results <= 2:
                    print('\n')
                    print(f'\nHi ! NCBI returned >{line[9:-9]}< which isn\'t \
enough for the subsequent analysis.I\'m sorry about this...\n Please try again!\n'.center(350, '-'))
                    
                    if yes_or_no('Would you like to enter a different protein family?') == 'y':
                        protein = enter_protein()
                    
                elif 2 < prot_tax_results < 1000:
                    print('\n')
                    print(f'\nThis program has found >{line[9:-9]}< results \
within your specified taxa sub-group. We can now continue with the analysis...\n'.center(350, '-'))
                    break
                        
                elif prot_tax_results >= 1000:
                    print('\n')
                    print('\nThis is a lot of sequences to get from NCBI and \
so this may take a while...\n'.center(200, '-')) 
                          
                    if yes_or_no('Would you still like to continue?') == 'y':
                        print('\n')
                        print('\nMake yourself a cup of tea!\n'.center(100, '-'))
                        break
                    
                    elif yes_or_no('Would you like to enter a different protein family?')=='y':
                        protein = enter_protein()
        
                else:
                    print("Something is going wrong, please try again!")
                                        
        return protein, taxonomic_subset

# =============================================================================
# Getting rid of partial sequences
# =============================================================================
def prune_seq(protein, organism):
    
    print('\nRetrieving sequences...\n'.center(120, '-'))
    
    sequence_cmd =f'esearch -db protein -query "{protein}[PROT] AND {organism}[ORGN]"\
    |efetch -format fasta>file1'
        
    subprocess.call(sequence_cmd, shell=True)
    
    print('\nPruning sequences...\n'.center(120, '-'))
    
    fasta_header_cmd = "grep '>' file1 > headers.txt"
    subprocess.call(fasta_header_cmd, shell=True)
    
    partial_cmd = 'grep -v "partial" headers.txt | cut --complement -c 1 > headers_wo_partial.txt'
    subprocess.call(partial_cmd, shell=True)
        
    pullseq_cmd = './pullseq -i file1 -n headers_wo_partial.txt > file_wo_partial'
    subprocess.call(pullseq_cmd, shell=True)
        
    count_cmd = 'seqcount file_wo_partial -outfile seqcount.out'
    subprocess.call(count_cmd, shell=True)
        
    cat_cmd = 'cat seqcount.out'
    count = subprocess.check_output(cat_cmd, shell=True)
    
    if int(count) > 2:
        
        return(int(count))
            
    else:
        print('\nI am so sorry it seems that after removing partial sequences returned by NCBI \
we don\'t have enough sequences to do a multiple alignment...\n'.center(350, '-'))
            
        if yes_or_no('Would you like to start again?') == 'n':
            print('\nexiting...\n'.center(250, '~'))
            
        else:
            print('\nStarting again...\n'.center(250, '~'))
            user_input()
            
# =============================================================================
# Number of species        
# =============================================================================
def species():
    
    with open('headers_wo_partial.txt') as myfile:
        lines = [line.rstrip() for line in myfile]
    
    species_subset = []
    
    
    for line in lines:
        species_subset.append(line.split('[')[1][:-1])
    
    species_number = len(set(species_subset))
    species_set = set(species_subset)
    
    return lines, species_number, species_set
    
            
# =============================================================================
# program order
# =============================================================================

def user_input():
    try:
        while True:
#            subprocess.call('clear', shell=True)
            prot=enter_protein()
            prot2, orgn=enter_organism(prot)
            
            try:
            
                number = prune_seq(prot2, orgn)
            except:
                print('prune_seq')
                
            try:
                headers, num_spec, set_spec = species()
            except:
                print('species')
        
            print(f'\nAccording to your input you want to look at the \
conservation of >{prot2}s< in >{num_spec}< different species... \nPlease be aware \
with many different species in your analysis you may not get very satisfactory\
 conservation results...\n'.center(450, '-'))
            
            if yes_or_no(f'Would you like to see the {num_spec} species?') == 'y':
                for i in set_spec:
                    print(i)
            
            if yes_or_no(f'With this in mind, would you like to continue with \
>{num_spec}< species?') == 'y':

                print('\nYay! Let\'s continue with the analysis...\n'.center(200, '~'))
                
                return headers, prot2, orgn
            
            else:
                if yes_or_no('Would you like to start again?') == 'n':
                    print('\nexiting...\n'.center(250, '~'))
                    break
                else:
                    print('\nStarting again...\n'.center(250, '~'))
                    
    except:
        print('USER INPUT FUNCTION')

# =============================================================================
# CHOOSING SEQUENCES FROM dict for PATMATMOTIF        
# =============================================================================
def choose_from_dict(question, dictionary):
    while True:
        num_choice=input("\n###\n"+question+"\n###\n")
        
        valid_integers = list(dictionary.keys())
        
        if int(num_choice) in valid_integers:
            break
        else:
            print('\nYou have not enetered a valid number\nTry again!\n'.center(300,'-'))
    return num_choice, dictionary[int(num_choice)]





    
