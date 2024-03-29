#!/usr/bin/python

__author__      = "Pierre Lesturgie"
__version__     = "0.3"
__email__       = "pierrelesturgie@outlook.fr"
__date__        = "2023-01-16"

import argparse
from pysam import VariantFile
from tqdm import tqdm

ap = argparse.ArgumentParser(description="Filter vcf for heterozygosity and write a filtered VCF")
ap.add_argument("-i", "--input", required=True,
                help="input vcf file must be compressed with bgzip and indexed with tabix (or uncompressed)")
ap.add_argument("-H", "--het", required=False,
                help="threshold for heterozygote frequency (default is 0.8)",type=float, default=0.8)
ap.add_argument("-O", "--output", required=False,
                help="output prefix",type=str, default="./output")
ap.add_argument("-c", "--chr", required=False,
                help="chromosome",type=str,default=None)
ap.add_argument("-d", "--write_discarded", required=False,
                help="write discarded SNPs (default is False)",action='store_true', default=False)

#print("WARNING: This version uses a stric threshold for heterozygosity and allows to write discarded SNPs")
#print(" >>> use v.0.2 for '<=' ")

args = vars(ap.parse_args())

# Reading input vcf file and opening output

vcf_in = VariantFile(str(args['input']))
vcf_out = VariantFile(str(args['output']) + ".vcf", 'w', header=vcf_in.header)

if args['write_discarded'] is True : 
    file = open(str(args['output']) + "_discarded_het.bed", 'w')

for rec in tqdm(vcf_in.fetch(args['chr'])):
    alleles = []
    het = 0
    for s in rec.samples.values() :
        if s['GT'][0] == None :
            alleles.append(s['GT'][0])
        elif s['GT'][0] <= 1 and s['GT'][1] <= 1:
            alleles.append(s['GT'][0]) 
            alleles.append(s['GT'][1])
            if s['GT']==(0,1) or s['GT']==(1,0):
                het = het + 1
    #print(het)
    if len(alleles) > 0:
        alleles = [int(0) if v is None else v for v in alleles]
        #print(het)
        if (het / (len(alleles) / 2)) <= args['het']:
            vcf_out.write(rec)
        elif args['write_discarded'] is True:
            file.write(rec.chrom + "\t" + str((rec.pos + 1)))
            file.write('\n') 
if args['write_discarded'] is True : 
    file.close()       
    
vcf_out.close()