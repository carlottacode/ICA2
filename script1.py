#!/usr/bin/python3
import os, sys, subprocess
import numpy as np
import pandas as pd

subprocess.call('echo "welcome to the ICA2"', shell=True)
subprocess.call("esearch -db protein -query 'pyruvate dehydrogenase[PROT]'|\
efilter -organism 'ascomycete fungi'|efetch -format fasta>file1", shell=True)
