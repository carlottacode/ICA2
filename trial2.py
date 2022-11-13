#!/usr/bin/python3

import os, sys, subprocess
import pandas as pd
import numpy as np

from functions2 import enter_protein
from functions2 import enter_organism
from functions2 import yes_or_no
from functions2 import species
from functions2 import user_input

#try:
# =============================================================================
#   #Error trapping protein family user input...
# =============================================================================
    #try:
        #prot=enter_protein()
    #except:
        #print('PROTEIN FUNCTION')
# =============================================================================
#   #Error trapping taxonomic subset user input...
# =============================================================================
    #try:
        #prot2, orgn=enter_organism(prot)
        #print(prot2)
        #print(orgn)
    #except:
        #print('ORGANISM FUNCTION')

#except:
    #print('WHOLE THING')

#try:
	#species(prot2, orgn)

#except:
	#print("SPECIES FUNCTION")

try:
	a, b=user_input()
	print(a+' '+b)
except:
	print('CONSERVATION')
