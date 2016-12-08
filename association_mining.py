"""

Building association rules into k-mer finding pipeline


export PATH=/Users/cbazodi/anaconda3/bin:$PATH

"""
import sys, os
import pandas as pd


def R_arules(DF, wd, num, m):
  """ Run arules.R script to make association rules and save as df_rules.csv """
  R=('R --vanilla --slave --args '+wd+' '+DF+' '+num+' '+m+' < /Users/cbazodi/Desktop/Combined_stress/ARuleMining/Arules.R')
  os.system(R)


def rule_test(df, row, r):
  """ Returns a Series for the rule given with
  1 = all kmers in rule present, 0 = at least 1 kmer in rule missing """

  rule_len = len(r)
  count = 0
  for ri in r:
    if row[ri] == 1:
      count += 1
    
  if rule_len == count:
    return 1
  else:
    return 0



def make_features(DF, wd, rules):
  """ Reads in the rule.csv and calls rule_test() to build new association rule features """
  
  # Make dictionary of all rules - lhs is/are the antecedent(s) and rhs is the consequence
  rule_dict = []
  for l in open(rules, 'r').readlines():
    if 'lhs' in l and 'rhs'  in l:
      pass
    else:
      line = l.replace('}', '').replace('{', '').strip().split('"')
      lhs = line[1].strip().split(",")
      rhs = line[5]
      r_kmers = lhs + [rhs]
      rule_dict += [r_kmers]
  print("Number of rules being converted to features: %i" % (len(rule_dict)))


  # 
  df = pd.read_csv(DF, sep='\t', index_col = 0)
  row = pd.Series(df.index.values)

  for r in rule_dict:
    feat_name = ':'.join(r)
    df[feat_name] = df.apply(lambda row: rule_test(df, row, r), axis=1)

  df2 = pd.concat([df['Class'], df.filter(like=':')], axis = 1)

  save_name1 = DF + "_plusARs"
  save_name2 = DF + "_onlyARs"
  
  df.to_csv(save_name1, sep="\t", float_format='i')
  df2.to_csv(save_name2, sep="\t", float_format='i')




if __name__ == "__main__":  

  if len(sys.argv) <= 1:
    print(__doc__)
    exit()

  n = 10
  m = 'lift'
  
  for i in range (1,len(sys.argv),2):
    if sys.argv[i] == "-df":            # data frame file (Col 1 = ID, Col2 = Class, Col3-... = Features)
      DF = sys.argv[i+1]
    if sys.argv[i] == "-wd":            # working directory
      wd = sys.argv[i+1]
    if sys.argv[i] == "-n":             # number of rules to mine
      num = sys.argv[i+1]
    if sys.argv[i] == "-m":             # measure for selecting best rules ('support', 'confidence','lift')
      m = sys.argv[i+1]

  R_arules(DF, wd, num, m)

  rules = wd+DF+"_rules.csv"

  rule_dict = make_features(DF, wd, rules)


  
