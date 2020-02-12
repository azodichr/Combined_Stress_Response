"""make_runcc_OverlapFilter.py
"""


import os, sys


for i in range (1,len(sys.argv),2):

  if sys.argv[i] == '-dir':    #Path to directory with files
    path = sys.argv[i+1]


other_path = "/mnt/home/azodichr/01_CombinedStress/Niu_MgCO2/03_Kmers/01_FS_tree"
script = "/mnt/home/azodichr/GitHub/Combined_Stress_Response/FilterDF_Overlap.py"
fun = "cat"
export = "export PATH=/mnt/home/azodichr/miniconda3/bin:$PATH"

out = open("runcc_CategAttributes.txt", 'w')

d = os.listdir(path)

for i in d:
  if i.startswith("._"):
    pass

  elif "imp" in i or "results" in i or "CNSorDHS" in i:
    pass

  elif "CNS" in i:
    df1 = os.path.join(path, i)
    name = i[:-10]
    name2 = name+"_NiFST_DHS"
    df2 = os.path.join(path, name2)
    df3 = os.path.join(other_path, name)
    save = name + "_Categorical_PA_CNS_DHS"
    out.write("%s; python %s -f %s -df %s -df2 %s -df3 %s -save %s\n" % (export, script, fun, df1, df2, df3, save))
