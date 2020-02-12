#make_mapping_runcc.py

import os, sys
for i in range (1,len(sys.argv),2):

  if sys.argv[i] == '-dir':    #Path to directory with files
    path = sys.argv[i+1]
  if sys.argv[i] == '-save':
    save = sys.argv[i+1]

out_name = "runcc_RF_" + save + ".txt"
out = open(out_name, 'w')

export = "PATH=/mnt/home/azodichr/miniconda3/bin:$PATH"
script = "/mnt/home/azodichr/GitHub/MotifDiscovery/RF_scikit.py"


for f in os.listdir(path):
  if f.startswith(".") or "imp" in f or "results" in f or "runcc" in f:
    pass
  elif "df" in f:
    name = f + "_" + save
    out.write("export %s; python %s -df %s -save %s\n" % (export, script, path+f, name))
