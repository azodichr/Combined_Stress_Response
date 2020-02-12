"""make_runcc_Mapped2GFF.py
"""


import os, sys

save = "skip"
script = "/mnt/home/azodichr/GitHub/Utilities/gff_manager.py"
fun = "mapped2gff"
key = "/mnt/research/ShiuLab/15_Dnaseseq/coordinates/COL.coord.1000promoter_all"
source = "RF_Kmer_Pipeline"
Type = "mapped_kmer"

for i in range (1,len(sys.argv),2):

  if sys.argv[i] == '-dir':    #Path to directory with files
    path = sys.argv[i+1]
  elif sys.argv[i] == '-source':    #Path to directory with files
    source = sys.argv[i+1]


out = open("runcc_Mapped2GFF.txt", 'w')

d = os.listdir(path)

for i in d:
  if i.startswith("._"):
    pass

  elif "out.pvalue" in i:
    x = os.path.join(path, i)
    out.write("python %s -f %s -i %s -key %s -source %s -type %s\n" % (script, fun, x, key, source, Type))
