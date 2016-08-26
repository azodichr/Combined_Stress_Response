import pandas as pd
import sys, os

alig = sys.argv[1]
key = sys.argv[2]

al = pd.read_csv(alig, sep = "\t", header=0, index_col = 0)
k = pd.read_csv(key, sep="\t", header = 0, index_col = 0, usecols = [0,1])

print(k.head(5))
print(al.head(3))

full = pd.merge(al, k, left_index=True, right_index=True, how = 'left')

print(full.head(2))

full.to_csv("agilent_021169_20160208_TAIR10_mod.txt", sep = "\t")