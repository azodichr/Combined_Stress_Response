#make_mapping_runcc.py

import os, sys

for i in range (1,len(sys.argv),2):

  if sys.argv[i] == '-dir':    #Path to directory with files
    path = sys.argv[i+1]
  if sys.argv[i] == '-include':
    inc = sys.argv[i+1]

out_name = "runcc_make_hybrid_df_FET.txt"
out = open(out_name, 'w')

script = "/mnt/home/azodichr/GitHub/Combined_Stress_Response/combined_df.py"
dap = "/mnt/home/azodichr/01_CombinedStress/peaks_parsed_ampRemoved.csv.mapped.regulator.matrix.altsplicRemoved.txt"
cis = "/mnt/home/azodichr/01_CombinedStress/01_ProcessingCISbp/CISBP_Hits2NonOverlappingPromoters_AT.matrix.txt"


for f in os.listdir(path):
  if f.startswith(".") or "imp" in f or "results" in f:
    pass
  elif inc in f:
    name = "hybridDF_" + f
    out.write("python %s -dap %s -cis %s -kmer %s -o %s\n" % (script, dap, cis, path+f, name))
