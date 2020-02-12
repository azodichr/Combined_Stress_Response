#make_mapping_runcc.py

import os, sys

for i in range (1,len(sys.argv),2):

  if sys.argv[i] == '-dir':    #Path to directory with files
    path = sys.argv[i+1]
  if sys.argv[i] == '-method':
    method = sys.argv[i+1]
  if sys.argv[i] == '-n':
    n = sys.argv[i+1]

out_name = "runcc_FeatSel_" + method + ".txt"
out = open(out_name, 'w')

export = "PATH=/mnt/home/azodichr/miniconda3/bin:$PATH"
script = "/mnt/home/azodichr/GitHub/MotifDiscovery/Feature_Selection_sklearn.py"


for f in os.listdir(path):
  if f.startswith(".") or "imp" in f or "results" in f:
    pass
  elif f.endswith("_df.txt"):
    out.write("export %s; python %s -f %s -df %s -n %s\n" % (export, script, method, path+f, n))
