"""make_runcc_OverlapFilter.py
"""


import os, sys

f_type = "gff"
cutoff = "na"


for i in range (1,len(sys.argv),2):

  if sys.argv[i] == '-dir':    #Path to directory with files
    path = sys.argv[i+1]
  if sys.argv[i] == '-save':    #Suffix for output files
    save = sys.argv[i+1]
  if sys.argv[i] == '-overlap':
    overlap = sys.argv[i+1]
  if sys.argv[i] == '-ftype':
    f_type = sys.argv[i+1]
  if sys.argv[i] == '-cutoff':
    cutoff = sys.argv[i+1]


script = "/mnt/home/azodichr/GitHub/Combined_Stress_Response/FilterDF_Overlap.py"
fun = "filter"
export = "export PATH=/mnt/home/azodichr/miniconda3/bin:$PATH"

out_save_name = "runcc_OverlapFilter_" + save + ".txt"
out = open(out_save_name, 'w')

d = os.listdir(path)

for i in d:
  if i.startswith("._"):
    pass

  elif "imp" in i or "results" in i or "FST_CNS" in i:
    pass

  elif "df" in i:
    if F_TYPE == "gff":
      x = os.path.join(path, i)
      out.write("%s; python %s -f %s -df %s -overlap %s -save %s\n" % (export, script, fun, x, overlap, save))
    
    elif F_TYPE == "bed":
      x = os.path.join(path, i)
      out.write("%s; python %s -f %s -df %s -overlap %s -save %s -ftype %s -cutoff %s\n" % (export, script, fun, x, overlap, save, f_type, cutoff))
