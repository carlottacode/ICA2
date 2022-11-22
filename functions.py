#!/usr/bin/python3
import os, sys, subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pprint
import textwrap


output_path = './prog_run/'


# =============================================================================
# These functions control the formatting of all the print statements for the 
# user. They simply take the string to be printed and wrap the text and 
# allow the messages to be differentiated using symbols.
# =============================================================================
def print_info(msg):
    print('\n')
    print(''.center(70, '-'))
    print(textwrap.fill(msg, 70))
    print(''.center(70, '-'))
    
def print_process(msg):
    print('\n')
    print(''.center(70, '/'))
    print(textwrap.fill(msg, 70))
    print(''.center(70, '\\'))
    
def print_checkpoint(msg):
    print('\n')
    print(''.center(70, '~'))
    print(textwrap.fill(msg, 70))
    print(''.center(70, '~'))
    
def print_graph(msg):
    print('\n')
    print(''.center(70, '.'))
    print(textwrap.fill(msg, 70))
    print(''.center(70, '.'))
    
# =============================================================================
# Create output file
# =============================================================================
def log_counter_init():
    exec_list = []
    if not os.path.exists('./log_file.txt'):
        
        os.mknod('./log_file.txt')
        os.chmod('./log_file.txt', 744)
    
    with open('log_file.txt', 'r+') as myfile:
        log_string = [line for line in myfile]
        print(log_string)
        if log_string == []:
            myfile.write('Execution: 1\n')
            log_counter = 1
            return log_counter
        else:
            for i in log_string:
                if 'Execution:' in i:
                    exec_list.append(i)
                    
            log_counter = int(exec_list[-1].split(' ')[1]) + 1
            myfile.write('Execution: '+str(log_counter)+"\n")
            return log_counter

# =============================================================================
# write to output file
# =============================================================================

def write_to_log(thing):
    with open('log_file.txt', 'a') as myfile:
        myfile.write(thing)

# =============================================================================
# This function allows you to ask a yes or no question which either returns a 
# counter 'y' or 'n'. I use this in other functions when deciding to interact 
# with the user. I have usually used this to exit a while True: loop
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
            print_info('I\'m sorry I don\'t understand... Please enter (y\\n) in the accepted format')

# =============================================================================
# This function allows the user choose from multiple options at once using a 
# dictionary and integers as keys. This function should only allow the user to 
# continue when they have entered a valid key.        
# =============================================================================
def check_dict(question, dictionary):
    while True:
        num_choice=input("\n###\n"+question+"\n###\n")
        
        valid_integers = list(dictionary.keys())
        try:
            if int(num_choice) in valid_integers:
                return num_choice, dictionary[int(num_choice)]
            else:
                print_info('You have not enetered a valid number. Try again!')
                
        except ValueError:
            print_info('You have not input a valid integer.')
                        
# =============================================================================
# This function allows the user to input a number, either a float or an integer
# in a specific range of numbers. This function should only allow the user to 
# continue when they have entered a number or float within in the range specified.
# =============================================================================
def check_num(max_len, min_len, question):
      while True:
        num_choice=input("\n###\n"+question+"\n###\n")
        try:
                num = int(num_choice)
                if (min_len < int(num_choice) < max_len) == True:
                      return num
                else:
                      print_info(f'The sequence length you have entered is not in the correct range should be between >{max_len}< and >{min_len}<')            
        except ValueError:
                try:
                    num = float(num_choice)
                    if (min_len < num < max_len) == True:
                          return num
                    else:
                         print_info(f'The sequence length you have entered is not in the correct range should be between >{max_len}< and >{min_len}<')
                except ValueError:
                    print_info('You have not entered a number, please try again...')
        
# =============================================================================
# This function makes sure that the user enters a string when prompted for 
# protein family and taxonomic subset.
# =============================================================================

def check_string(question):
    while True:
        query = input("\n###\n"+question+"\n###\n").lower()
        if query.replace(' ', '').isalpha():
            return query
        else:
            print_info("You have entered some weird characters - enter characters A-Z only")

# =============================================================================
# This function asks the user to enter a protein family of interest. It will 
# not allow the user to continue unless the protein family returns sequences
# from NCBI. Therefore hopefully catching all user input errors.
# =============================================================================
def enter_protein():
        while True:
            
                protein_family=check_string('What protein family are you interested in?')
                #protein_family = 'pyruvate dehydrogenase'
                
                if protein_family[-1] == 's':
                  protein_family = protein_family[:-1]
                  
                #Checking that the protein family returns results from NCBI if not then 
                #there may be spelling mistakes 
                esearch_cmd = f'esearch -db IPG -query "{protein_family}[PROT]" | grep "Count"'
                
                line =  subprocess.check_output(esearch_cmd, shell=True).decode("utf-8")
                protein_results = int(line[9:-9])
                
                if protein_results == 0:
                    print_info(f">{protein_family[0].upper()+protein_family[1:]}< didn\'t return any results from NCBI...This means you might have made a typo... Please try again!")
                              
                elif 1<= protein_results <= 2:
                    print_info(f'Hi ! NCBI returned >{protein_results}< which isn\'t enough for the subsequent analysis. I\'m sorry about this... Please try again!')
                        
                elif protein_results > 2 :
                    print_info(f"Luckily for you NCBI has returned >{line[9:-9]}< sequence results for your protein family (>{protein_family}<)!")
                    break
                else:
                    print_info("Something isn't working please try again...")
        return protein_family


# =============================================================================
# This functions asks the user to enter a taxonomic subset. It doesn't allow 
# the user to continue with the analysis until a combination of protein family 
# and taxonomic subset returns more than 2 sequences from NCBI. There are points
# at which if the number of sequences returned is not satisfactory that you are
# able to enter a different protein family if you wish, after which 
# you can enter another taxonomic sub-group.
# =============================================================================
def enter_organism(protein):
        while True:
                taxonomic_subset=check_string("What taxonomic subset do you wish to investigate?")
                
                efilter_cmd = f'esearch -db IPG -query "{protein}[PROT] AND {taxonomic_subset}[ORGN]"|grep "Count"'
                
                line =  subprocess.check_output(efilter_cmd, shell=True).decode("utf-8")
                prot_tax_results = int(line[9:-9])
                
                if prot_tax_results == 0:
                    print_info(f"Hello, your combination of >{taxonomic_subset}< and >{protein}< didn\'t return any results from NCBI...This means you might have made a typo... Please try entering the taxonomic subgroup again!")
                    
                    if yes_or_no('Would you like to enter a different protein family?') == 'y':
                        protein = enter_protein()
                              
                elif 1 <= prot_tax_results <= 2:
                    print_info(f'Hi ! NCBI returned >{prot_tax_results}< which isn\'t enough for the subsequent analysis.I\'m sorry about this... Please try again!')
                    
                    if yes_or_no('Would you like to enter a different protein family?') == 'y':
                        protein = enter_protein()
                    
                elif 2 < prot_tax_results < 1000:
                    print_info(f'This program has found >{prot_tax_results}< results within your specified taxa sub-group. We can now continue with the analysis...')
                    break
                        
                elif prot_tax_results >= 1000:
                    print_info(f'Processing >{prot_tax_results}< may take a while...')
                          
                    if yes_or_no('Would you still like to continue?') == 'y':
                        print_info('Make yourself a cup of tea!')
                        break
                    
                    elif yes_or_no('Would you like to enter a different protein family?')=='y':
                        protein = enter_protein()
        
                else:
                    print_info("Something is going wrong, please try again!")
                                        
        return protein, taxonomic_subset

# =============================================================================
# This function retrieves the fasta sequences for chosen the protein and 
# taxonomic sub-group. The fasta headers are put in another file. Which is 
# searched for partial sequences. A file is created with headers of only the complete
# sequences. A list of these headers is returned.
# =============================================================================
def prune_seq(protein, organism, path=output_path):
    
    print_process('Retrieving sequences...')
    sequence_cmd =f'esearch -db protein -query "{protein}[PROT] AND {organism}[ORGN]"|efetch -format fasta>{path}file1'
    subprocess.call(sequence_cmd, shell=True)
    
    print_process('Removing partial sequences...')
    fasta_header_cmd = f"grep '>' {path}file1 > {path}headers.txt"
    subprocess.call(fasta_header_cmd, shell=True)
    
    partial_cmd = f'grep -i -v "partial" {path}headers.txt | cut --complement -c 1 > {path}headers_wo_partial.txt'
    subprocess.call(partial_cmd, shell=True)
    
    with open(f'{output_path}headers_wo_partial.txt') as myfile:
        lines = [line.rstrip() for line in myfile]
        
    return lines
            
# =============================================================================
# This function creates a set of the species returned for the chosen protein 
# and taxonomic subgroup using another NCBI search.
# =============================================================================
def species(prot, orgn):
    species_cmd = f"esearch -db IPG -query '{prot}[PROT] AND {orgn}[ORGN]'|esummary|xtract -pattern DocumentSummary -element Organism"
    species_list = subprocess.check_output(species_cmd, shell=True).decode("utf-8")
    species_subset=species_list.split('\n')
    species_number = len(set(species_subset))
    species_set = set(species_subset)
    
    return species_number, species_set

# =============================================================================
# program order
# =============================================================================

def user_input(path=output_path):
    try:
        while True:
            subprocess.call('clear', shell=True)
            prot=enter_protein()
            prot2, orgn=enter_organism(prot)
            
            headers = prune_seq(prot2, orgn)
            
            num_spec, set_spec = species(prot2, orgn)
            
            # If there are enough complete sequences to continue 
            # with the analysis pullseq is used to create a clustalo input file which 
            # contains only the complete sequences returned by NCBI. If there are not enough
            # sequences the user is prompted whether they would like to start again which
            # would allow them to input a new protein group and subsequently a new taxonomic
            # sub-group.
            if len(headers) > 2 :
                
                pullseq_cmd = f'./pullseq -i {path}file1 -n {path}headers_wo_partial.txt > {path}file_wo_partial'
                subprocess.call(pullseq_cmd, shell=True)
                
                print_info(f'According to your input you want to look at the conservation of >{prot2}s< in >{num_spec}< different species... Please be aware with many different species in your analysis you may not get very satisfactory\conservation results...')
                
                if yes_or_no(f'Would you like to see the >{num_spec}< species?') == 'y':
                    for i in set_spec:
                        print(i)
                
                if yes_or_no(f'With this in mind, would you like to continue with >{num_spec}< species?') == 'y':
                    print_checkpoint('Yay! Let\'s continue with the analysis...')
                    return headers, prot2, orgn
                    break 
                
                else:
                    print_info('I am so sorry it seems that after removing partial sequences returned by NCBI we don\'t have enough sequences to do a multiple alignment...')
                    
                    if yes_or_no('Would you like to start again?') == 'n':
                        print_checkpoint('exiting...')
                        break
                    else:
                        print_checkpoint('Starting again...')
                        
            else:
                if yes_or_no('Would you like to start again?') == 'n':
                    print_checkpoint('exiting...')
                    break
                else:
                    print_checkpoint('Starting again...')
                    
    except:
        print('THE USER INPUT FUNCTION FAILED')
                
# =============================================================================
# This function provides the user with the option to use a subset of sequences 
# in the multiple sequence alignment instead of all sequences. 
# The user is presented with the sequence dataframe and a histogram to get an
# idea of how many sequences there are of a particular length. 
# The user can then either remove short or long sequences or choose their own
# range entirely. 
# The user can only progress if their range returns more than 2 sequences for
# the multiple sequence alignment.
# =============================================================================
def define_subset(df):
    print_graph('Plotting a histogram of the sequence lengths')
    
    plt.hist(df['Sequence Length'], color=('lightgreen'), ec='darkgreen')
    plt.xlabel('Sequence Length')
    plt.ylabel('Number of Sequences')
    plt.title('Sequence Lengths for your Search Query')
    plt.savefig(output_path+'seq_histogram.pdf', format='pdf')
    
    cmd_view_hist = f'gs {output_path}seq_histogram.pdf'
    subprocess.call(cmd_view_hist, shell=True)
    
    
    my_max = df['Sequence Length'].loc[df['Sequence Length'].idxmax()]      # Maximum in column
    my_min = df['Sequence Length'].loc[df['Sequence Length'].idxmin()]      # Maximum in column
    
    num_seq = len(df.index.tolist())
    
    if yes_or_no('Would you like to specify a subset?') == 'n':
        print_info(f'Continuing with all >{num_seq}< sequences...')
        seq = df.index.tolist()
    
    else: 
        
        while True:
            subset_options = {1:'Specify a subset by removing some of the shorter sequences?', 
                              2:'Specify a subset by removing some of the longer sequences?',
                              3:'Specify a subset by defining your own range?'}
            
            
            pprint.pprint(subset_options)
            key, dict_input=check_dict('How would you like to create your subset?', subset_options)
            
            print('\n')
            pprint.pprint(subset_options)
            
            if key == '1':
                num = check_num(my_max, my_min, 'Please enter a minimum sequence length.')
            
                subset=df[df['Sequence Length']>num]
                num_seq = len(subset.index.tolist())
                
                seq = subset.index.tolist()
                
                if num_seq > 2:
                    if yes_or_no(f'You have chosen {num_seq} sequences. Would you like to choose a different subset?') == 'n':
                        print_checkpoint(f'Continuing with this subset of >{num_seq}< sequences!')
                        break
                
                else:
                    print_info(f'I\'m really sorry but your subset only contains >{num_seq}< sequences and you need at least three to to a multiple sequence alignemnt')
            
                
            if key == '2':
                num = check_num(my_max, my_min, 'Please enter a maximum sequence length.')
            
                subset=df[df['Sequence Length']<num]
                num_seq = len(subset.index.tolist())
                seq = subset.index.tolist()
                
                
                if num_seq > 2:
                    if yes_or_no(f'You have chosen {num_seq} sequences. Would you like to choose a different subset?') == 'n':
                        print_checkpoint(f'Continuing with this subset of >{num_seq}< sequences!')
                        break
                else:
                    print_info(f'I\'m really sorry but your subset only contains >{num_seq}< sequences and you need at least three to to a multiple sequence alignment')
                
            if key == '3':
                min_num = check_num(my_max, my_min, 'Please enter your minimum sequence length.')
                max_num = check_num(my_max, my_min, 'Please enter your maximum sequence length.')
                
                subset=df[df['Sequence Length'].between(min_num, max_num)]
                num_seq = len(subset.index.tolist())
                seq = subset.index.tolist()
                
                
                if num_seq > 2:
                    if yes_or_no(f'You have chosen {num_seq} sequences. Would you like to choose a different subset?') == 'n':
                        
                        print_checkpoint(f'Continuing with this subset of >{num_seq}< sequences!')
                        break
                else:
                    print_info(f'I\'m really sorry but your subset only contains >{num_seq}< sequences and you need at least three to to a multiple sequence alignemnt')
            
    return num_seq, seq
            
# =============================================================================
# This function carries out a multiple sequence alignment using Clustal Omega. 
# Calls a function allowing the user to specify a subset of sequences to 
# align and then uses pullseq to create an input file with only the sequences 
# chosen by the user. Clustalo is then run creating specific output files 
# which are used for further analysis.
# =============================================================================
def clustalo(protein, organism, path=output_path):
    if yes_or_no(f'Would you like to align the sequences you retrived from NCBI for >{protein}< in >{organism}<?') == 'y':
        
        with open(f'{output_path}file_wo_partial') as myfile:
            seq_lines= [line for line in myfile]
        
        seq =[]  
        for i in seq_lines:
            if '>' in i:
                seq.append(i.replace('\n','>'))
            else:
                seq.append(i)
            
        seq=''.join(seq)
        seq = seq.replace('X', '').replace('\n', '').split('>')[1:]
        seq_df = pd.DataFrame(seq[1::2], index=seq[::2], columns=['Sequence'])
        seq_df['Sequence Length'] = seq_df['Sequence'].str.len()
        seq_df.sort_values(['Sequence Length'], ascending=False, inplace=True)
        
        print(seq_df)
        
        #Using my function define_subset to return sequences the user wants to look at by length.
        num_seq, subset_seq = define_subset(seq_df)
                    
        
        
        if yes_or_no('Would you like to set a maximum number of sequences to take from this subset?')=='y':
            max_seqs=check_num(num_seq, 2, 'How many of these sequences would you like to do a multiple sequence alignment on?')
        else:
            max_seqs=num_seq
        
        
        seq_file      = []
        for i in subset_seq[0:max_seqs]:
            i = i+'\n'
            seq_file.append(i)               
        
        
        with open(f'{path}msa_in', 'w') as myfile:
            myfile.writelines(seq_file)
        
        pullseq_cmd = f'./pullseq -i {path}file1 -n {path}msa_in>{path}clustalo_in'  
        subprocess.call(pullseq_cmd, shell = True)
        cmd = f"clustalo --force --full --threads 16 --outfmt=phy --maxnumseq={max_seqs} --guidetree-out={path}clustalo.tree -i {path}clustalo_in -o {path}clustalo.phy"
        subprocess.call(cmd, shell=True)
        
        
# =============================================================================
# HMOMENT
# =============================================================================
def hmoment(seq_dict, path=output_path):
    print(seq_dict)
    plt.figure(figsize=(30, 10))
    counter = 0
    while True:
        
        if yes_or_no('Would you like to add a protein sequence to your hydropathy plot?')=='y':
            counter += 1
            pprint.pprint(seq_dict)
            
            key, dict_input=check_dict('To create a hydropathy plot enter the number of the sequence', seq_dict)        
        
            temp_file_cmd=f'echo {dict_input}>{path}hp_temp_file'
            pullseq_cmd = f'./pullseq -i {path}file1 -n {path}hp_temp_file > {path}hp_temp_seqs'
        
            subprocess.call(temp_file_cmd, shell=True)
            subprocess.call(pullseq_cmd, shell=True)
            
            hp_cmd = f'hmoment {path}hp_temp_seqs -outfile {path}hmom_data'
            subprocess.call(hp_cmd, shell=True)
            
            with open(f'{path}hmom_data') as myfile:
                to_plot = myfile.read()
                
            hmom_plot = to_plot.split('HMOMENT')[1:]
            
            leg_label='HMOMENT'+hmom_plot[0].split('\n')[0]    
            data = [y.split('\t\t') for y in hmom_plot[0].split('\n')[4:-1]]
            df = pd.DataFrame(data, columns=hmom_plot[0].split('\n')[3].split('\t'))
            df = df.apply(pd.to_numeric)
            plt.plot(df['Position'], df['uH'], label=leg_label)
            
            plt.xlabel('Position')
            plt.title('Hydropathy Plots of Protein Sequences')
            plt.ylabel('uH')
            plt.legend(loc=1)
            plt.savefig(f'{path}hydropathy_plot.pdf', format='pdf')
            
            view_hp_cmd = f'gs {path}hydropathy_plot.pdf'
            subprocess.call(view_hp_cmd, shell=True)
            write_to_log(f'Hydropathy plot for {dict_input}: hydropathy_plot.pdf\n')
        else:
            if counter<1:
                print_info('We can\'t make a plot without any sequences!')
            
            else:
                break
        
