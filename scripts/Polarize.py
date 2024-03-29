#!/usr/bin/python

__author__      = "Pierre Lesturgie"
__version__     = "0.1"
__email__       = "pierrelesturgie@outlook.fr"
__date__        = "2023-04-12"


import argparse
from pysam import VariantFile, tabix_compress, tabix_index
import itertools
from tqdm import tqdm
from subprocess import run

ap = argparse.ArgumentParser(description="polarize a VCF given an outgroup")

ap.add_argument("-i", "--input", required=True,
                help="input vcf file file")
ap.add_argument("-o", "--output", required=True,
                help="output prefix for vcf")
ap.add_argument("-O", "--outgroup", required=True,
   help="Name of the outgroup individual -- MUST be included in the VCF :-) ", type=str)
ap.add_argument("-D", "--write_derived_site", required=False,
   help="Write a table with the number of derived alleles per individual per site", 
                action='store_true', default=False)
ap.add_argument("-d", "--write_derived_ind", required=False,
   help="Write a table with the total number of derived alleles per individual ** This is 5x longer", 
                action='store_true', default=False)
ap.add_argument("-t", "--tempDir", required=False,
                help="temporary directory",default="tmp", type=str)
ap.add_argument("-m", "--includeMiss", required=False,
                help="do not filter out missing data for vcf and .der file", 
                action='store_true', default=False)
ap.add_argument("-w", "--include_outgroup", required=False,
                help="Include OG individual in the VCF", 
                action='store_true', default=False)

args = vars(ap.parse_args())

run("mkdir " + str(args['tempDir']),shell=True)

temp = open(str(args['tempDir']) + "/temp.vcf", 'w')
temp.write("##fileformat=VCFv4.2\n")
temp.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
temp.close()

run("zgrep contig " + str(args['input']) + " > "  + str(args['tempDir']) + "/temp",shell=True)
run("cat " + str(args['tempDir']) + "/temp.vcf " + str(args['tempDir']) + "/temp > " + str(args['tempDir']) + "/temp_head.vcf",shell=True)
#run("rm temp.vcf temp",shell=True)

f = VariantFile(str(args['input']), 'r')

vcf = open(str(args['tempDir']) + "/temp_gen.vcf", 'w')
SAMPLES=list(f.header.samples)
OG=str(args['outgroup'])
index_outgroup=SAMPLES.index(OG)

if args['include_outgroup'] is False :
    SAMPLES.pop(index_outgroup)

coln=[["#CHROM", "POS", "ID","REF", "ALT","QUAL","FILTER","INFO","FORMAT"]]
coln.append(SAMPLES)
coln = [k for i in coln for k in i]
vcf.write('\t'.join(str(i) for i in coln))
vcf.write('\n')

if args['write_derived_site'] is True or args['write_derived_ind'] is True: 
    calcul_DER = True
    if args['write_derived_site'] is True : 
        DER = open(args['output'] + ".der", 'w')
        g = [["CHROM", "POS"]]
        g.append(SAMPLES)
        g = list(itertools.chain(*g))
        DER.write('\t'.join(str(i) for i in g))
        DER.write('\n')
    if args['write_derived_ind'] is True : 
        derind = open(str(args['tempDir']) + "/count.derind", 'w')
        derind.write('\t'.join(str(i) for i in [0]*len(SAMPLES)))
        derind.write('\n')
        derind.close()
        tempspl = open(str(args['tempDir']) + "/sample.derind", 'w')
        tempspl.write('\t'.join(str(i) for i in SAMPLES))
        tempspl.write('\n')
        tempspl.close()
else : 
    calcul_DER = False
    


for rec in tqdm(f.fetch()):
    MD_FLG = False
    alleles = []
    for s in rec.samples.values() :
        alleles.append(s.alleles)
    if args['include_outgroup'] is False :
        og_alleles=alleles.pop(index_outgroup)
    else : 
        og_alleles = alleles[index_outgroup]
    if og_alleles[0] == og_alleles[1] and og_alleles[0] is not None: 
        REF=og_alleles[0]
        ALT=rec.ref if rec.ref != REF else rec.alts
        if len(list(itertools.chain(*ALT))) == 1 and len(list(itertools.chain(REF))) == 1: 
            genotypes = [rec.contig,rec.pos,".",REF,ALT[0],rec.qual,".",".","GT"]
            der = [rec.contig,rec.pos]
            for g in alleles: 
                if g[0] == REF and g[1] == REF:
                    genotypes.append("0/0")
                    if calcul_DER is True:
                        der.append("0")
                elif g[0] != REF and g[1] != REF:
                    if g[0] is None:
                        genotypes.append("./.")
                        if calcul_DER is True:
                            der.append(".")
                            MD_FLG=True
                    else:
                        genotypes.append("1/1")
                        if calcul_DER is True:
                            der.append("2")
                else:
                    genotypes.append("0/1")
                    if calcul_DER is True:
                        der.append("1")
            if args['tempDir'] is True or MD_FLG is not True:
                vcf.write('\t'.join(str(i) for i in genotypes))
                vcf.write('\n')
            if calcul_DER is True:
                if args['tempDir'] is True or MD_FLG is not True: 
                    if args['write_derived_site'] is True : 
                        DER.write('\t'.join(str(i) for i in der))
                        DER.write('\n')
                    if args['write_derived_ind'] is True : 
                        derind = open(str(args['tempDir']) + "/count.derind", 'r')
                        keep = []
                        for line in derind: 
                            k = [int(i) for i in line.split('\t')]
                        derind.close()
                        derind = open(str(args['tempDir']) + "/count.derind", 'w')
                        d = [int(i) for i in der[2:]]
                        ret = [x + y for x, y in zip(k, d)]
                        derind.write('\t'.join(str(i) for i in ret))
                        derind.write('\n')
                        derind.close()
                
                
vcf.close()

if args['write_derived_site'] is True : 
    DER.close()

if args['write_derived_ind'] is True : 
    run("cat " + str(args['tempDir']) + "/sample.derind " + str(args['tempDir']) + "/count.derind > " + str(args['output']) + ".derind",shell=True)
    #run("rm sample.derind count.derind",shell=True)

run("cat " + str(args['tempDir']) + "/temp_head.vcf " + str(args['tempDir']) + "/temp_gen.vcf > " + str(args['output']) + ".vcf",shell=True)
#run("rm temp_gen.vcf temp_head.vcf",shell=True)

tabix_compress(str(args['output']) + ".vcf",  str(args['output']) + ".vcf.gz",force=True)
tabix_index(str(args['output']) + ".vcf.gz",preset='vcf',force=True)
run("rm " + str(args['output']) + ".vcf",shell=True)
                 
run("rm -r " + str(args['tempDir']), shell = True)