""" Extract kmers out of hit lists """

import os
import pandas as pd


d1 = "/mnt/home/azodichr/01_CombinedStress/Niu_MgCO2/03_Kmers/01_FS_tree"
d2 = "/mnt/home/azodichr/01_CombinedStress/Perata_HAnox/03_Kmers/01_FS_tree"
d3 = "/mnt/home/azodichr/01_CombinedStress/Prasch_HD/03_Kmers/01_FS_tree"
d4 = "/mnt/home/azodichr/01_CombinedStress/Niu_MgCO2/03_Kmers/02_FS_Fishers"
d5 = "/mnt/home/azodichr/01_CombinedStress/Perata_HAnox/03_Kmers/02_FS_Fishers"
d6 = "/mnt/home/azodichr/01_CombinedStress/Prasch_HD/03_Kmers/02_FS_Fishers"

directories = [d1, d2, d3, d4, d5, d6]
kmers = []

for d in directories:
  d_i = os.listdir(d)
  
  for i in d_i:
    if i.startswith("._"):
      pass

    elif "_imp.txt" in i:
      x = os.path.join(d, i)
      f = pd.read_csv(x, sep='\t', index_col = 0)
      kmer_list = list(f.index)
      kmers.extend(kmer_list)

print("Total number of kmers (including replicates between clusters): %i" % len(kmers))

kmers_noreps = list(set(kmers))

print("Number of kmers with replicates removed: %i" % len(kmers_noreps))

out = open("kmers_list.txt", 'w')

for k in kmers_noreps:
  out.write("%s\n" % k)