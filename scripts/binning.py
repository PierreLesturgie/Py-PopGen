#!/usr/bin/python

__author__      = "Pierre Lesturgie"
__version__     = "2.0"
__email__       = "pierrelesturgie@outlook.fr"
__date__        = "2023-04-15"

import pandas as pd
import argparse
from pysam import VariantFile
from tqdm import tqdm

# Construct the argument parser
ap = argparse.ArgumentParser(description="Perform a binning from a VCF and write the binned VCF")

# Add the arguments to the parser
ap.add_argument("-i", "--input", required=True,
                help="input vcf file must be compressed with bgzip and indexed with tabix (or uncompressed)")
ap.add_argument("-o", "--output", required=True,
                help="output name for vcf")
ap.add_argument("-b", "--bin", required=True,
                help="bin size")
ap.add_argument("-r", "--region", required=False,
                help="region size",default=None)
args = vars(ap.parse_args())



def binning(vcf_in, vcf_out,bin_length,reg_length=None) : 
    vcf_in = VariantFile(str(vcf_in))
    vcf_out = VariantFile(str(vcf_out), 'w', header=vcf_in.header)
    bin_length = int(bin_length)
    
    bound_inf = 0
    chro = list((vcf_in.header.contigs))[0]
    
    if reg_length is not None:
        reg_length = int(reg_length) 
        bound_supp = bound_inf + reg_length

    for rec in tqdm(vcf_in.fetch()):
        if rec.chrom == chro:
            if reg_length is None:
                if bound_inf == 0 or rec.pos > bound_inf + bin_length:
                    vcf_out.write(rec)
                    bound_inf = rec.pos
            else : 
                if rec.pos >= bound_inf and rec.pos < bound_supp: 
                    vcf_out.write(rec)
                if rec.pos >= bound_supp + bin_length: 
                    bound_inf = bound_supp + bin_length
                    bound_supp = bound_inf + reg_length

        else: 
            chro = rec.chrom
            bound_inf = 0
            if reg_length is not None:
                bound_supp = bound_inf + reg_length 

    vcf_out.close()

    
binning(vcf_in=args['input'], vcf_out=args['output'],
        bin_length=args['bin'],reg_length=args['region'])