"""
Combined DAP-Seq, CIS-bp, and kmer dataframes for machine learning

"""

import pandas as pd
import numpy as np
import sys

for i in range (1,len(sys.argv),2):
  if sys.argv[i] == "-dap":
    DAP = sys.argv[i+1]
  if sys.argv[i] == "-cis":
    CIS = sys.argv[i+1]
  if sys.argv[i] == "-kmer":
    KMER = sys.argv[i+1]
  if sys.argv[i] == "-o":
    save = sys.argv[i+1]

# Set up dataframe
kmer = pd.read_csv(KMER, sep='\t', index_col = 0)
cis = pd.read_csv(CIS, sep='\t', index_col = 0)
dap = pd.read_csv(DAP, sep='\t', index_col = 0)

# Join all three on genes from kmer dataframe (has the class variable)
two = kmer.join(dap)
three = two.join(cis)

# Save final combined dataframe
save_name = save + "_df.txt"
three.to_csv(save_name, sep="\t", na_rep = "0")