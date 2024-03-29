#!/usr/bin/python

__author__      = "Pierre Lesturgie"
__version__     = "0.3"
__email__       = "pierrelesturgie@outlook.fr"
__date__        = "2022-07-12"

from pysam import VariantFile
import pandas as pd
import numpy as np
import itertools
from tqdm import tqdm
import argparse
import subprocess

# Construct the argument parser
ap = argparse.ArgumentParser(description='Compute Unphased Extended Haplotype using r2')

# Add the arguments to the parser
ap.add_argument("-i", "--input", required=True,
   help="input vcf file if compressed must be with bgzip and indexed with tabix")
ap.add_argument("-b", "--bin", required=True,
   help="bin bounds : <from> <to> <by>", nargs='+', type=int)
ap.add_argument("-w", "--window", required=True,
   help="window size", type=int)
ap.add_argument("-s", "--bootstrap", required=False,
   help="number of bootstrap", type=int, default=None)
ap.add_argument("-r", "--runs", required=True,
   help="number of runs per window", type=int)
ap.add_argument("-f", "--buffer", required=False,
   help="buffer size", type=int, default = 1000)
ap.add_argument("-m", "--maf", required=True,
   help="maf filter", type=float)
ap.add_argument("-c", "--chr", required=False,
    help="chromosome",type=str,default=None)
ap.add_argument("-t", "--region", required=False,
   help="region bounds : <from> <to>", nargs='+', type=int, default=None)
args = vars(ap.parse_args())

#print(args["window"])
#print(args["input"])



def extended_geno(pos_df,gen_df,window,bootstrap,buff,bins):
    np.seterr(invalid='ignore')
    for l in tqdm(range(0,pos_df.col.to_list()[-1],window)):
        #print("<<< Computing extended haplotype for " + lines[l].split()[0] + " to " + lines[l].split()[1] + " window")
        if l != 0:
            win_max = l
            candidates = pd.DataFrame({'col':pos_df.index[pos_df["col"].between(win_min,win_max)].tolist()})
            r2_data_frame = []
            if len(candidates.col) > 0:
                for j in bootstrap:
                    SNP = list(itertools.chain(*(candidates.sample().values.tolist())))
                    pairs = []
                    for i in bins :
                        sup, inf = int(pos_df.loc[SNP[0]]) + i, int(pos_df.loc[SNP[0]]) - i
                        bound_sup = pd.DataFrame({'col':pos_df.index[pos_df["col"].between(sup - buff,sup + buff)].tolist()})
                        bound_inf = pd.DataFrame({'col':pos_df.index[pos_df["col"].between(inf - buff,inf + buff)].tolist()})
                        if len(bound_sup) > 0 : 
                            bound_sup = list(itertools.chain(*(bound_sup.sample().values.tolist())))
                            pairs.append([[SNP[0]],bound_sup,[i]])
                        if len(bound_inf) > 0 : 
                            bound_inf = list(itertools.chain(*(bound_inf.sample().values.tolist())))
                            pairs.append([[SNP[0]],bound_inf,[-i]])
                    r2 = []
                    for k in range(len(pairs)):
                        temp = list(itertools.chain(*pairs[k]))
                        gen_temp = gen_df.iloc[temp[0:2]]
                        delta = int(pos_df.loc[temp[1]]) - int(pos_df.loc[temp[0]])
                        g0 = list(itertools.chain(*(list(itertools.chain(*gen_temp.iloc[0])))))
                        g1 = list(itertools.chain(*(list(itertools.chain(*gen_temp.iloc[1])))))
                        interval = [str(win_min) + "_" + str(win_max)]
                        interval.append(int(pos_df.loc[SNP[0]]))
                        interval.append(delta)
                        interval.append(temp[2])
                        interval.append(np.corrcoef(g0, g1)[0,1]**2)
                        interval.append(int(j))
                        r2.append(interval)
                    r2_data_frame.append(pd.DataFrame(r2))
                r2_data_frame = pd.concat(r2_data_frame)
                r2_data_frame.to_csv("temp_dir_EH/" + str(win_min) + "_" + str(win_max) + "_Extended_Haplotype.ld",index=False,header=False,sep=' ')
        win_min = l


# Reading input VCF file
vcf_in = VariantFile(str(args['input']))

print("<<< formating genotype and filtering maf >>>")

if args['region'] is None : 
    START, STOP = None, None
else : 
    START, STOP = args['region'][0], args['region'][1]

maf, genotypes, positions = [], [], []
for rec in tqdm(vcf_in.fetch(args['chr'], START, STOP)):
    alleles = []
    geno_temp = []
    for s in rec.samples.values() :
        alleles.append(s['GT'][0])
        alleles.append(s['GT'][1])
        if s['GT'] == (0,0) : 
            geno_temp.append([0])
        elif s['GT'] == (1,1) :
            geno_temp.append([1])
        else : 
            geno_temp.append([2])
    maf_temp = sum(alleles) / len(alleles)  
    maf.append(maf_temp)
    if maf_temp >= float(args['maf']) and maf_temp <= (1 - float(args['maf'])) :
        genotypes.append(geno_temp)
        positions.append(rec.pos)

pos_df, gen_df = pd.DataFrame({'col':positions}), pd.DataFrame({'col':genotypes})

print("<<< computing extended haplotype r2 value >>>")

subprocess.run(["mkdir", "temp_dir_EH"])

extended_geno(pos_df=pos_df,gen_df=gen_df,window=args['window'],
              bootstrap=list(range(0,args["runs"],1)),buff=args["buffer"],
              bins=list(range(args["bin"][0],args["bin"][1],args["bin"][2])))

if args["bootstrap"] is not None :
    print("<<< computing bootstraps >>>")
    r2 = []
    r2_data_frame = []
    for i in range(args["bootstrap"]):
        SNP = list(itertools.chain(*(candidates.sample(2).values.tolist())))
        gen_temp = gen_df.iloc[SNP]
        g0 = list(itertools.chain(*(list(itertools.chain(*gen_temp.iloc[0])))))
        g1 = list(itertools.chain(*(list(itertools.chain(*gen_temp.iloc[1])))))
        interval = [np.corrcoef(g0, g1)[0,1]**2]
        interval.append(int(i))
        r2.append(interval)
        r2_data_frame.append(pd.DataFrame(r2))
    r2_data_frame = pd.concat(r2_data_frame)
    r2_data_frame.to_csv("bootstrap.ld",index=False,header=False,sep=' ')


names=['window','SNP_Pos','Delta','bins','r2',"bootstrap"]
subprocess.run(["cat temp_dir_EH/*.ld > temp_dir_EH/Extended_Haplotype"],shell=True)
pd.DataFrame(names).transpose().to_csv("temp_dir_EH/names.names",header=False,sep=' ',index=False)
subprocess.run(["cat temp_dir_EH/names.names temp_dir_EH/Extended_Haplotype > ./Extended_Haplotype.ld"],shell=True)
subprocess.run(["rm -r temp_dir_EH"],shell=True)