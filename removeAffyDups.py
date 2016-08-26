#script_removeAffyDups.py
#Remove lines from affy ATH1 probe code if the probe hits multiple loci.

import os, sys
import pandas as pd

f = sys.argv[1]

df = pd.read_csv(f, sep = '\t', header=0, index_col = 0)

df = df[df.locus.str.contains(";") == False]

df.to_csv("affy_ATH1_array_elements-2010-12-20_noDups.txt", sep = "\t")