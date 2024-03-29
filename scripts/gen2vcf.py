#!/usr/bin/python

__author__      = "Pierre Lesturgie"
__version__     = "0.1"
__email__       = "pierrelesturgie@outlook.fr"
__date__        = "2023-02-24"

import pandas as pd
import argparse
from tqdm import tqdm

ap = argparse.ArgumentParser(description="Converts FSC2 .GEN file to vcf")

ap.add_argument("-i", "--input", required=True,
                help="input .gen file")
ap.add_argument("-o", "--output", required=True,
                help="output name for vcf")
ap.add_argument("-C", "--chr_info", required=True,
   help="number of chromosome and length simulated : <nchr> <len>", nargs='+', type=int)

args = vars(ap.parse_args())


f = open(str(args['input']), 'r')
vcf = open(str(args['output']), 'w')
vcf.write("##fileformat=VCFv4.2\n")
vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
nchr=args['chr_info'][0]
lenchr=args['chr_info'][1]

for i in range(nchr):
    #print(i)
    vcf.write("##contig=<ID=" + str(i+1) + ",length=" + str(lenchr) + ",assembly=SIMULATION.fa>\n")



for line in tqdm(f):
    if line[0] == "C" : 
        tempo = line.split()
        colnames = [tempo.pop(0),tempo.pop(0),"ID",tempo.pop(0),tempo.pop(0),"QUAL","FILTER","INFO","FORMAT"]
        colnames[0], colnames[1], colnames[3], colnames[4] = "#CHROM", "POS", "REF", "ALT"
        while tempo: 
            temp = [tempo.pop(0),tempo.pop(0)]
            colnames.append(temp[0])
            #ind_name = temp[0].split("_")[0] + "_" + temp[0].split("_")[1]
        vcf.write('\t'.join(str(i) for i in colnames))
        vcf.write('\n')
    else: 
        row = line.split()
        if len(row)>0: 
            genotypes = [row.pop(0),row.pop(0),".",row.pop(0),row.pop(0),".",".",".","GT"]
            row = [int(v) for v in row]
            while row: 
                temp = [row.pop(0),row.pop(0)]
                if temp[0]!=temp[1]:
                    temp2 = "0/1"
                else : 
                    if temp[0] == 1:
                        temp2 = "1/1"
                    else : 
                        temp2 = "0/0"
                genotypes.append(temp2)
            vcf.write('\t'.join(str(i) for i in genotypes))
            vcf.write('\n')
vcf.close()