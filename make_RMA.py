import sys, os
import pandas as pd


for i in range (1,len(sys.argv),2):
  if sys.argv[i] == '-dir':             # pwm sequence file (TAMO)
    d = sys.argv[i+1]
  if sys.argv[i] == '-o':             # pwm sequence file (TAMO)
    o = sys.argv[i+1]

names = []
n = 0
for f in os.listdir(d):
  if f.startswith("."):
    pass
  else:
    n += 1
    name = f.strip().split("_")[1]
    names.append(name)
    path = d + f
    
    if n == 1:
      matrix = pd.read_csv(path, sep = "\t", skiprows = 2, header=None)

    else:
      add = pd.read_csv(path, sep = "\t", skiprows = 2, header=None)
      matrix = pd.merge(matrix, add, how = 'inner', on = 0)

matrix = matrix.set_index(0)
matrix.columns = names
print(matrix.head(5))
print(matrix.shape)

matrix_na = matrix.dropna(axis=0, how='any')
print(matrix_na.shape)

out_name = o + "_rma.txt"
matrix_na.to_csv(out_name, sep = "\t")
