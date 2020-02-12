"""make_runcc_CombineFiltered.py
"""

import os, sys


for i in range (1,len(sys.argv),2):

  if sys.argv[i] == '-dir':    #Path to directory with files
    path = sys.argv[i+1]


script = "/mnt/home/azodichr/GitHub/Combined_Stress_Response/FilterDF_Overlap.py"
fun = "combine"
export = "export PATH=/mnt/home/azodichr/miniconda3/bin:$PATH"

out = open("runcc_CombineFilter.txt", 'w')

d = os.listdir(path)

files = []

for i in d:
  if i.startswith("._"):
    pass

  elif "imp" in i or "results" in i:
    pass

  elif "CNS" in i or "DHS" in i:
    files.append(i)

for f in files:
  if "CNS" in f:
    df1 = f
    search = df1[:-4]
    for x in files:
      if search in x and "DHS" in x:
        df2 = x
    save = search + "_CNSorDHS"
    print(df1,df2)

    df1x = os.path.join(path, df1)
    df2x = os.path.join(path, df2)

    out.write("%s; python %s -f %s -df %s -df2 %s -save %s\n" % (export, script, fun, df1x, df2x, save))
