#!/usr/bin/python

__author__      = "Pierre Lesturgie"
__version__     = "0.2"
__email__       = "pierrelesturgie@outlook.fr"
__date__        = "2023-04-05"

import argparse
from pysam import VariantFile
import pandas as pd
import numpy as np
import itertools
from tqdm import tqdm

# Construct the argument parser
ap = argparse.ArgumentParser()

# Add the arguments to the parser
ap.add_argument("-i", "--input", required=True,
   help="input vcf")
ap.add_argument("-o", "--output", required=True,
   help="output name")
ap.add_argument("-b", "--bin", required=True,
   help="bin size")
ap.add_argument("-f", "--buffer", required=True,
   help="buffer size")
ap.add_argument("-m", "--maf", required=True,
   help="minor allele frequency")
args = vars(ap.parse_args())

# Reading input VCF file
vcf = VariantFile(str(args['input']))

# Extracting Genotypes, Positions
# Computing and filtering maf

print(" >>> reading genotypes, positions and filtering maf...")
maf, genotypes, positions = [], [], []
none = 0
for rec in tqdm(vcf.fetch()):
    alleles = []
    geno_temp = []
    for s in rec.samples.values() :
        alleles.append(s['GT'][0])
        alleles.append(s['GT'][1])
        if s['GT'] == (0,0) : 
            geno_temp.append([0])
        elif s['GT'] == (None,None) : 
            geno_temp.append([None])
        elif s['GT'] == (1,1) :
            geno_temp.append([1])
        else : 
            geno_temp.append([2])
    if len([value for value in alleles if value is None]) > 0:
        none += 1
    else : 
        maf_temp = sum(alleles) / len(alleles)  
        maf.append(maf_temp)
        if maf_temp >= float(args['maf']) and maf_temp <= (1 - float(args['maf'])) :
            genotypes.append(geno_temp)
            positions.append(rec.pos)

pos_df, gen_df = pd.DataFrame({'col':positions}), pd.DataFrame({'col':genotypes})
print("Discarded", none, "sites with missing data")
print("...Done. <<< ")


print(" >>> Finding pairs of SNPs for LD using a bin size of", int(args['bin']), "bp and a buffer of", int(args['buffer']), "bp...")
pairs = []
# Finding the pairs on which to compute LD
for i in tqdm(positions): 
    mini = i + int(args['bin']) - int(args['buffer']) 
    maxi = i + int(args['bin']) + int(args['buffer']) 

    bound_sup = pd.DataFrame({'col':pos_df.index[pos_df["col"].between(mini,maxi)].tolist()})

    if len(bound_sup) > 0 : 
        bound_sup = list(itertools.chain(*(bound_sup.sample().values.tolist())))
        bound_inf = pos_df.index[pos_df['col']==i].tolist()
        pairs.append([bound_inf,bound_sup])
        #print((i/positions[-1])*100,"%")
print("...Done. <<< ")


print(" >>> Computing genotypic r2...")
# Computing genotypic r2
r2 = []
np.seterr(invalid='ignore')
for i in tqdm(range(len(pairs))): 
    temp = list(itertools.chain(*pairs[i]))
    gen_temp = gen_df.iloc[temp]
    pos_temp = pos_df.iloc[temp]
    g0 = list(itertools.chain(*(list(itertools.chain(*gen_temp.iloc[0])))))
    g1 = list(itertools.chain(*(list(itertools.chain(*gen_temp.iloc[1])))))
    interval = list(itertools.chain(*[pos_temp.iloc[0].tolist(),pos_temp.iloc[1].tolist()]))
    interval.append(np.corrcoef(g0, g1)[0,1]**2)
    r2.append(interval)
    #print((i/(len(pairs)-1))*100,"%")
print("...Done. <<< ")

r2_df = pd.DataFrame(r2, columns=['Pos1','Pos2','r2'])
pd.DataFrame.to_csv(r2_df,str(args['output']) + ".csv",sep=';',index=False)

print(r2_df.iloc[:, 2::3].describe())