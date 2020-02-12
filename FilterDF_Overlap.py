"""Functions with different methods to incorporate DAP-Seq data into machine learning pipelines.
-f:
    filter_by_overlap_gff:    
    DAP_enrich:     

"""


from collections import defaultdict
import pandas as pd
import sys, os
from datetime import datetime


class DAP_Seq:

  def filter_by_overlap_gff(self, OVERLAP, DF, SAVE):
    """Given a presence/absence dataframe this script removes presence hits for gene-kmer pairs that exist, but 
    don't overlap with a DAP region (or DNaseI region, basically whatever is in your overlap file!"""

    #Load dataframe for machine learning and save feature (i.e. kmer) names 
    if isinstance(DF, str):                               #If loading df as a file
      df = pd.read_csv(DF, sep='\t', index_col = False)   
    else:                                                 #If df is already in pandas format from previous script
      df = DF

    genes_in_matrix = df.iloc[:,0]
    wanted_kmers = list(df.columns.values)[2:]

    poss = df.loc[df['Class'] == 1]
    pos_pre_filt = sum(poss[wanted_kmers].sum(axis=1))
    negs = df.loc[df['Class'] == 0]
    neg_pre_filt = sum(negs[wanted_kmers].sum(axis=1))

    print(pos_pre_filt, neg_pre_filt)
    print("Genes in matrix: " + str(len(genes_in_matrix)) + " \nKmers in matrix: %" + str(len(wanted_kmers)))
    
    ##Make filter_dict which contains all the kmers in your dataframe and the genes in which the kmer overlaps a DAP region/motif site in the gene promoter.
    # {kmer: ['AT2G37260', 'AT2G37260', 'AT2G37260', 'AT2G37260', etc]}
    filter_dict = defaultdict(list)
    count_noOverlap = 0
    for l in open(os.path.join(OVERLAP),'r'):
      line = l.strip().split("\t")
      kmer = line[2]
      gene = line[1]
      if kmer in wanted_kmers:
      # gene in genes_in_matrix:
        attributes = line[8].split(";")
        overl = attributes[2]
        if overl == "overlap=NA":
          pass
          count_noOverlap = count_noOverlap + 1
        else:
          filter_dict[kmer].append(gene)
      
      else:
        pass

    print("Number of kmers with overlaps: " + str(len(filter_dict))) 

    count_progress_1 = 0
    for kmer_col in wanted_kmers:
      
      if count_progress_1 % 10 == 0:
        print ("Completed filtering presence/absense matrix for " + str(count_progress_1) + " kmers")
      else:
        pass
      count_progress_1 = count_progress_1 + 1

      for idx in range(len(df)):
        g = df.iloc[idx][0]
        if g in filter_dict[kmer_col]:
          #print(g + " matches " + kmer_col)
          df[kmer_col].iloc[idx] = 1
        else:
          df[kmer_col].iloc[idx] = 0

    poss_post = df.loc[df['Class'] == 1]
    pos_post_filt = sum(poss_post[wanted_kmers].sum(axis=1))
    negs_post = df.loc[df['Class'] == 0]
    neg_post_filt = sum(negs_post[wanted_kmers].sum(axis=1))

    save_name = DF + "_" + SAVE
    df.to_csv(save_name, sep="\t", index = False, float_format='%.0f')
    open("/mnt/home/azodichr/01_CombinedStress/Filtering_Summary.txt",'a').write('%s\t%i\t%i\t%i\t%i\n' % (DF, pos_pre_filt, neg_pre_filt, pos_post_filt, neg_post_filt))
    

  def filter_by_overlap_bed(self, OVERLAP, CUTOFF, DF, SAVE):
    """Given a presence/absence dataframe this script removes presence hits for gene-kmer pairs that exist, but 
    don't overlap with a DAP region (or DNaseI region, basically whatever is in your overlap file!"""

    #Load dataframe for machine learning and save feature (i.e. kmer) names 
    if isinstance(DF, str):                               #If loading df as a file
      df = pd.read_csv(DF, sep='\t', index_col = False)   
    else:                                                 #If df is already in pandas format from previous script
      df = DF

    genes_in_matrix = df.iloc[:,0]
    wanted_kmers = list(df.columns.values)[2:]

    poss = df.loc[df['Class'] == 1]
    pos_pre_filt = sum(poss[wanted_kmers].sum(axis=1))
    negs = df.loc[df['Class'] == 0]
    neg_pre_filt = sum(negs[wanted_kmers].sum(axis=1))

    print(pos_pre_filt, neg_pre_filt)
    print("Genes in matrix: " + str(len(genes_in_matrix)) + " \nKmers in matrix: %" + str(len(wanted_kmers)))
    
    ##Make filter_dict which contains all the kmers in your dataframe and the genes in which the kmer overlaps a DAP region/motif site in the gene promoter.
    # {kmer: ['AT2G37260', 'AT2G37260', 'AT2G37260', 'AT2G37260', etc]}
    filter_dict = defaultdict(list)
    
    for l in open(os.path.join(OVERLAP),'r'):
      if l.startswith("Region"):
        pass
      else:
        line = l.strip().split("\t")    # line: [4_6327453_6327459_AT4G10150|6327453-6327459|AATATA, NoDiffSites, NoDiffSites]
        print(line)
        info = line[0].split("_")[3].split("|")    # info: [AT4G10150, 6327453-6327459, AATATA]
        kmer = info[2]
        gene = info[0]
        if kmer in wanted_kmers:
          conserved = line[1]
          if conserved == "NoDiffSites" or float(conserved) <= float(CUTOFF):
            filter_dict[kmer].append(gene)   


    print("Number of kmers with overlaps: " + str(len(filter_dict))) 

    count_progress_1 = 0
    for kmer_col in wanted_kmers:
      
      if count_progress_1 % 10 == 0:
        print ("Completed filtering presence/absense matrix for " + str(count_progress_1) + " kmers")
      else:
        pass
      count_progress_1 = count_progress_1 + 1

      for idx in range(len(df)):
        g = df.iloc[idx][0]
        if g in filter_dict[kmer_col]:
          #print(g + " matches " + kmer_col)
          df[kmer_col].iloc[idx] = 1
        else:
          df[kmer_col].iloc[idx] = 0

    poss_post = df.loc[df['Class'] == 1]
    pos_post_filt = sum(poss_post[wanted_kmers].sum(axis=1))
    negs_post = df.loc[df['Class'] == 0]
    neg_post_filt = sum(negs_post[wanted_kmers].sum(axis=1))

    save_name = DF + "_" + SAVE
    df.to_csv(save_name, sep="\t", index = False, float_format='%.0f')
    open("/mnt/home/azodichr/01_CombinedStress/Filtering_Summary.txt",'a').write('%s\t%i\t%i\t%i\t%i\n' % (DF, pos_pre_filt, neg_pre_filt, pos_post_filt, neg_post_filt))
    


  def combine_filtered(self, DF, DF2, SAVE):
    """Given filtered presence/absence dataframes this script adds them together so a 1 is given if the kmer overlaps
    with either criteria"""

    #Load dataframe for machine learning and save feature (i.e. kmer) names 
    df1 = pd.read_csv(DF, sep='\t', index_col = 0)   
    df2 = pd.read_csv(DF2, sep='\t', index_col = 0) 
    df_add = df1.add(df2, fill_value=0)
    df_add = df_add.replace(to_replace = 2, value = 1 )

    wanted_kmers = list(df1.columns.values)[2:]

    poss_post = df_add.loc[df_add['Class'] == 1]
    pos_post_filt = sum(poss_post[wanted_kmers].sum(axis=1))
    negs_post = df_add.loc[df_add['Class'] == 0]
    neg_post_filt = sum(negs_post[wanted_kmers].sum(axis=1))
    
    #print(pos_post_filt, neg_post_filt)
    df_add.to_csv(SAVE, sep="\t", float_format='%.0f')
    open("/mnt/home/azodichr/01_CombinedStress/Filtering_Summary.txt",'a').write('%s\t%s\t%i\t%i\n' % (DF, SAVE, pos_post_filt, neg_post_filt))

  def categorical_attributes(self, DF, DF2, DF3, SAVE):
    """Given dataframes for presense/absense, overlap with CNS, overlap with DHS = add them together
    0 = absent
    1 = present 
    2 = present + overlaps with CNS or DHS
    3 = present + overlaps with both CNS and DHS
    """

    #Load dataframe for machine learning and save feature (i.e. kmer) names 
    df1 = pd.read_csv(DF, sep='\t', index_col = 0)   
    df2 = pd.read_csv(DF2, sep='\t', index_col = 0) 
    df3 = pd.read_csv(DF3, sep='\t', index_col = 0)

    df_add = df1.add(df2, fill_value=0)
    df_add = df_add.add(df3, fill_value=0)

    df_add["Class"] = df_add["Class"].replace(to_replace = 3, value = 1 )

    df_add.to_csv(SAVE, sep="\t", float_format='%.0f')


  def DAP_enrich(self, OVERLAP, POS, NEG, SAVE, ALL):
    """Given parsed file of all genes that have a DAPpeak overlaping their x region (like promoter),
    this script makes a dataframe with all of the TFs that are enriched for hits in positive example genes
    """
    
    import numpy as np
    from scipy.stats import fisher_exact
    from math import sqrt

    pos_genes = [line.rstrip('\n').rstrip('\r') for line in open(POS)]
    neg_genes = [line.rstrip('\n').rstrip('\r') for line in open(NEG)]

    overlap = defaultdict(list)
    for l in open(OVERLAP,'r'):
      line = l.strip().split("\t")
      gene = line[0]
      tf = line[1].strip("Coords_").strip(".txt")
      overlap[tf].append(gene)

    TFs = overlap.keys()
    positive_present = {}.fromkeys(TFs, 0)     # Count occurence of each feature in positive examples
    negative_present = {}.fromkeys(TFs, 0)     # Count occurence of each feature in negative examples
    all_arrays = np.zeros([1,len(overlap.keys())+1])  # This fits the np df into the pd df - the plus 1 is for the Class! 
    

    #For each gene make an array with presence/absence of each feature in the overlap file
    for p in pos_genes:
      gene_array = np.array([1])
      for i in TFs:
        if p in overlap[i]:
          gene_array = np.append(gene_array, 1)
          positive_present[i] = positive_present[i]+1
        else:
          gene_array = np.append(gene_array, 0)
      all_arrays = np.vstack((all_arrays, gene_array))
    print("done with positives")
    
    for n in neg_genes:
      gene_array = np.array([0])
      for i in TFs:
        if n in overlap[i]:
          gene_array = np.append(gene_array, 1)
          negative_present[i] = negative_present[i]+1
        else:
          gene_array = np.append(gene_array, 0)
      all_arrays = np.vstack((all_arrays, gene_array))
    print("done with negatives")

    #Make the full dataframe
    columns = ['Class'] + list(overlap.keys())
    index = ['skip_this_line'] + pos_genes + neg_genes
    df= pd.DataFrame(all_arrays, index=index, columns=columns, dtype=int)  # Turn numpy into pandas DF
    df= df.drop("skip_this_line",0)
    
    if ALL == "true":
      save_name = SAVE + "_all_DAPtf_df.txt"
      df.to_csv(save_name, sep="\t", float_format='%.0f')


    # Calculate enrichement scores
    outFISH = open(SAVE+"_FETresults.txt",'w')
    outFISH.write('feature\tPosCount\tNegCount\tpvalue')
    count = 0
    
    enriched_kmers = {}
    for k in positive_present:
      try:
        count += 1
        TP = positive_present[k]            #Positive Examples with kmer present
        FP = negative_present[k]            #Negative Examples with kmer present
        TN = len(neg_genes)-negative_present[k]    #Negative Examples without kmer
        FN = len(pos_genes)-positive_present[k]    #Positive Examples without kmer

        oddsratio,pvalue = fisher_exact([[TP,FN],[FP,TN]],alternative='greater')
        outFISH.write('\n%s\t%d\t%d\t%.7f' % (k, (positive_present[k]),(negative_present[k]),pvalue))
        if pvalue <= PVAL:          # Remove unenriched features from dataframe
          enriched_kmers[k] = pvalue
        if pvalue > PVAL:
          df = df.drop(k, 1)
        if count%10000==0:
          print("Completed " + str(count) + " features")

      except ValueError:
        count += 1 
        outFISH.write('\n%s\t%d\t%d\t1.0' % (k, (positive_present[k]),(negative_present[k])))
    
    save_name2 = SAVE + "_enriched_" + str(PVAL) + "_DAPtf_df.txt"
    df.to_csv(save_name2, sep="\t", float_format='%.0f')





#-------------------------------------------------------------------------------

if __name__ == '__main__':
  DAP_Seq=DAP_Seq() # Print main function help if no imputs given
  ALL = "true"
  PVAL = 0.01
  F_TYPE = "gff"

  for i in range (1,len(sys.argv),2):

        if sys.argv[i] == '-overlap':    #Path to directory with DAP-Seq motifs mapped to peaks: /mnt/home/azodichr/DAPseq/00_raw_data/peaks/03_Map_PWM/01_StrictMapping/
          OVERLAP = sys.argv[i+1]
        if sys.argv[i] == '-ftype':     # Type of overlap file - optins = gff (default) or bed
          F_TYPE = sys.argv[i+1]
        if sys.argv[i] == '-cutoff':     # Type of overlap file - optins = gff (default) or bed
          CUTOFF = sys.argv[i+1]
        if sys.argv[i] == '-kmers':     #Path to directory with kmers mapped to promoter regions: /mnt/home/azodichr/01_DualStress_At/14_LogicClusters/06_MapKmers/01_tamo/01_LogicalCluster/
          KMERS = sys.argv[i+1]
        if sys.argv[i] == '-coords':    #Path to directory with kmers mapped to promoter regions: /mnt/home/azodichr/01_DualStress_At/14_LogicClusters/06_MapKmers/01_tamo/01_LogicalCluster/
          COORDS = sys.argv[i+1]
        if sys.argv[i] == '-df':        #Dataframe of pos-neg examples and kmer features with presence and absense of kmers.
          DF = sys.argv[i+1]
        if sys.argv[i] == '-df2':        #Dataframe of pos-neg examples and kmer features with presence and absense of kmers.
          DF2 = sys.argv[i+1]
        if sys.argv[i] == '-df3':        #Dataframe of pos-neg examples and kmer features with presence and absense of kmers.
          DF3 = sys.argv[i+1]
        if sys.argv[i] == '-pos':       #String to identify mapping files for positive examples
          POS = sys.argv[i+1]
        if sys.argv[i] == '-neg':       #String to identify mapping files for negative examples
          NEG = sys.argv[i+1]
        if sys.argv[i] == '-all':       #If TRUE, don't filter by enriched, keep all TFs
          ALL = sys.argv[i+1]
        if sys.argv[i] == '-save':       #String for save name
          SAVE = sys.argv[i+1]
        if sys.argv[i] == '-f':   
          F = sys.argv[i+1]


  if len(sys.argv) <= 1:
    print(__doc__)
    sys.exit()
  
  if F == "filter":
    if F_TYPE == "gff":
      if "" in [OVERLAP, DF, SAVE]:
        print("Need path to overlap file (from Johnny's gff_manager.py -f overlap) and to kmer p/a matrix to filter")
      DAP_Seq.filter_by_overlap_gff(OVERLAP, DF, SAVE)
    elif F_TYPE == "bed":
      if "" in [OVERLAP, CUTOFF, DF, SAVE]:
        print("Need path to overlap file (from Johnny's gff_manager.py -f overlap) and to kmer p/a matrix to filter, and a cutoff for what Nuc Diversity to count as conserved")
      DAP_Seq.filter_by_overlap_bed(OVERLAP, CUTOFF, DF, SAVE)

  if F == "combine":
    if "" in [DF, DF2, SAVE]:
      print("Need two kmer p/a matrixs to combine and a save suffix")
    DAP_Seq.combine_filtered(DF, DF2, SAVE)
  
  if F == "cat":
    if "" in [DF, DF2, DF3, SAVE]:
      print("Need two kmer p/a matrixs to combine and a save suffix")
    DAP_Seq.categorical_attributes(DF, DF2, DF3, SAVE)

  elif F == "DAP_enrich":
    if "" in [OVERLAP, POS, NEG, ALL, SAVE]:
      print("Need file with what genes have a DAP-peak overlapping a region (i.e. promoter), lists of positive and negative example genes, and a save name")
    DAP_Seq.DAP_enrich(OVERLAP, POS, NEG, SAVE, ALL)





