""" 

Converts .gff file to .bed file. Modified from Shin-Han's shius/codes/GFFUtil.py -f gff_to_bed

This version makes sure the kmer sequence is included in the description column

"""

import os, sys

for i in range (1,len(sys.argv),2):
        if sys.argv[i] == "-gff":
                gff = sys.argv[i+1]

inp = open(gff)
oup = open(gff+".bed","w")
inl = inp.readlines()
for i in inl:
        if i[0] != "#":
                L = i.strip().split("\t")
                C = L[0]        # chr
                cL= L[3]        # chr left
                cR= L[4]        # chr right
                #S = L[5]        # score
                try:
                        G = L[8].split(";")[0].split("=")[1] + "|" + L[2]
                        oup.write("%s\t%s\t%s\t%s\n" % (C,cL,cR,G))
                except IndexError:
                        print "Problem:",L[8]

print "Done!"