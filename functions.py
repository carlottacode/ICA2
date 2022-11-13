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
                protein_family=input("\nWhat protein family are you interested in?\n".center(123, '#')).lower()
                #protein_family = 'pyruvate dehydrogenase'
                
                if protein_family[-1] == 's':
                  protein_family = protein_family[:-1]
                  
                #Checking that the protein family returns results from NCBI if not then 
                #there may be spelling mistakes 
                esearch_cmd = f'esearch -db protein -query \
                    "{protein_family}[PROT]" | grep "Count"'
                
                line =  subprocess.check_output(esearch_cmd, shell=True).decode("utf-8")
                protein_results = int(line[9:-9])
                
                if protein_results == 0:
                    print(f"\nHi there! Your protein (>{protein_family}<) \
didn\'t return any results from NCBI...\nThis means you might have made a typo... \n\
Please try again!\n".center(123, '-'))
                              
                elif 1<= protein_results <= 2:
                    print(f'\nHi ! NCBI returned >{protein_results}< which \
isn\'t enough for the subsequent analysis. I\'m sorry about this... \n\
Please try again!\n'.center(123, '-'))
                        
                elif protein_results > 2 :
                        print(f"\nLuckily for you NCBI has returned \
>{line[9:-9]}< sequence results for your protein family (>{protein_family}<)!\n".center(123, '-'))
    
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
                taxonomic_subset=input("\nWhat taxonomic subset do you wish to investigate?\n".center(123, '#')).lower()
                #taxonomic_subset='ascomycete fungi'
                
                
                
                efilter_cmd = f'esearch -db protein -query "{protein}[PROT]"\
                    |efilter -organism "{taxonomic_subset}"|grep "Count"'
                
                line =  subprocess.check_output(efilter_cmd, shell=True).decode("utf-8")
                prot_tax_results = int(line[9:-9])
                
                if prot_tax_results == 0:
                    print(f"\nHello, your combination of >{taxonomic_subset}< \
and >{protein}< didn\'t return any results from NCBI... This means you might \
have made a typo... Please try entering the taxonomic subgroup again!\n".center(123, '-'))
                    
                    if yes_or_no('Would you like to enter a different protein family?') == 'y':
                        protein = enter_protein()
                              
                elif 1 <= prot_tax_results <= 2:
                    print(f'\nHi ! NCBI returned >{line[9:-9]}< which isn\'t \
enough for the subsequent analysis.I\'m sorry about this...\n Please try again!\n'.center(123, '-'))
                    
                    if yes_or_no('Would you like to enter a different protein family?') == 'y':
                        protein = enter_protein()
                    
                elif 2 < prot_tax_results < 1000:
                        print(f'\nThis program has found >{line[9:-9]}< results \
within your specified taxa sub-group. We can now continue with the analysis...\n'.center(123, '-'))
                        break
                        
                elif prot_tax_results >= 1000:
                    print('\nThis is a lot of sequences to get from NCBI and \
so this may take a while...\n'.center(123, '-')) 
                          
                    if yes_or_no('Would you still like to continue?') == 'y':
                        print('\nMake yourself a cup of tea!\n'.center(123, '-'))
                        break
                    
                    elif yes_or_no('Would you like to enter a different protein family?')=='y':
                        protein = enter_protein()
        
                else:
                    print("Something is going wrong, please try again!")
                                        
        return protein, taxonomic_subset




def species(protein, organism):
        try:
            sequence_cmd =f'esearch -db protein -query "{protein}[PROT]"|efilter -organism "{organism}"\
            |efetch -format fasta>file1'
                
            subprocess.call(sequence_cmd, shell=True)
              
    
            fasta_header_cmd = "grep '>' file1 > headers.txt"
            subprocess.call(fasta_header_cmd, shell=True)
        
        
            with open('headers.txt') as myfile:
                lines = [line.rstrip() for line in myfile]
            
            species_subset = []
            for line in lines:
                species_subset.append(line.split('[')[1][:-1])
                
            species_number = len(set(species_subset))
            
        except:
            print('SPECIES FUNCTION')
        
        return species_number

# =============================================================================
# PROGRAM 
# =============================================================================

def user_input():
    try:
        while True:
            prot=enter_protein()
            prot2, orgn=enter_organism(prot)    
            
            num1 = species(prot2, orgn)
        
            print(f'\nAccording to my records and the taxonomic subgroup you specified \
you want to look at the conservation of >{prot2}s< in >{num1}< different \
species... Please be aware with a large number of species in your analysis you \
may not get any cool or fun results...\n'.center(123, '-'))
            
            if yes_or_no(f'Baring this in mind, would you like to continue with \
>{num1}< species?') == 'y':

                print('...')
                return prot2, orgn
            else:
                if yes_or_no('Would you like to start again?') == 'n':
                    break
    except:
        print('USER INPUT FUNCTION')
        
