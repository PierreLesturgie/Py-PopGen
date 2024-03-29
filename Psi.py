#!/usr/bin/python

__author__      = "Pierre Lesturgie"
__version__     = "0.1"
__email__       = "pierrelesturgie@outlook.fr"
__date__        = "2023-04-15"

import argparse
from pysam import VariantFile, tabix_compress, tabix_index
import itertools
from tqdm import tqdm
from subprocess import run
import numpy

ap = argparse.ArgumentParser(description="compute Psi as in Peter & Slatkin (2013) and significancy")

ap.add_argument("-i", "--input", required=True,
                help="input .der file (** output from Polarize.py)")
ap.add_argument("-o", "--output", required=False,
                help="output prefix",default="out")
ap.add_argument("-r", "--resampling", required=False,
   help="Number of sfs resampling to perform [None]", default=None)
ap.add_argument("-t", "--tempDir", required=False,
                help="temporary directory",default="tmp", type=str)

args = vars(ap.parse_args())

run("mkdir " + str(args['tempDir']),shell=True)

run("cat " + str(args['input']) + " | wc -l > " + str(args['tempDir']) +  "/end.end", shell=True)

IN = open(str(args['input']),'r')

END = open(str(args['tempDir']) +  "/end.end",'r')

for line in END:
    end = int(line.splitlines()[0].strip())

print('\n' + " >>> Computing fixed shared alleles per pair of individuals" + "\n")
    
for line in tqdm(IN,total=end): 
    temp = line.splitlines()[0].split('\t')
    if temp[0] == "CHROM":
        keep_samples = temp[2:]
        p_samples = [str(i) for i in keep_samples]
        for i in range(len(p_samples)):
            for j in range(i,len(p_samples)):
                if j > i:
                    file_ext_temp = open(str(args['tempDir']) + "/" + p_samples[i] + "_" + p_samples[j] + ".temp", 'w')
                    file_ext_temp.write('\t'.join(str(i) for i in [0]*3))
                    file_ext_temp.close()
    else:
        p = [str(i) for i in temp[2:]]
        h = [1 for i in p if i != '.'] 
        if sum(h) == len(p): 
            k = [int(i) for i in p]
            for i in range(len(k)):
                for j in range(i,len(k)):
                    if j > i:
                        site, derived_i, derived_j = 1, 0, 0
                        if (k[i] == 2 or k[j] == 2) and k[i] != k[j]:
                            if k[i] != 0 and k[j] != 0 :
                                if k[i] == 2: 
                                    derived_i += 1   ### Het for j and fixed der for j
                                else : 
                                    derived_j += 1
                        file_ext_temp = open(str(args['tempDir']) + "/" + p_samples[i] + "_" + p_samples[j] + ".temp", 'r')
                        for line in file_ext_temp: 
                            keep = [int(i) for i in line.split('\t')]
                        file_ext_temp.close()
                        file_ext_temp = open(str(args['tempDir']) + "/" + p_samples[i] + "_" + p_samples[j] + ".temp", 'w')
                        d = [site, derived_i, derived_j]
                        ret = [x + y for x, y in zip(keep, d)]
                        file_ext_temp.write('\t'.join(str(i) for i in ret))
                        file_ext_temp.write('\n')
                        file_ext_temp.close()


print('\n' + " >>> Computing Psi (resampling = " + str(args['resampling']) + ") \n")

for i in tqdm(range(len(k))):
    for j in range(i,len(k)):
        if j > i:
            name = [[p_samples[i] + "_" + p_samples[j]]]
            file_ext_temp = open(str(args['tempDir']) + "/" + p_samples[i] + "_" + p_samples[j] + ".temp", 'r')
            for line in file_ext_temp: 
                keep = [int(i) for i in line.split('\t')]
            file_ext_temp.close()
            res = (keep[1] - keep[2])/keep[0]
            keep.append(res)
            if args['resampling'] is not None:
                boot=[]
                for b in range(int(args['resampling'])):
                    rng = numpy.random.default_rng()  
                    f = rng.binomial(n=1,p=0.5, size=keep[1] + keep[2])
                    boot.append((2*sum(f) - keep[1] - keep[2])/keep[0])
                a = [i for i in boot if i > res and res > 0]
                b = [i for i in boot if i < res and res < 0]
                prob = (len(a) + len(b)) / int(args['resampling'])
            else: 
                prob = "NaN"
            keep.append(prob)
            name.append(keep)
            res_fin = list(itertools.chain(*name))
            file_ext_temp = open(str(args['tempDir']) + "/" + p_samples[i] + "_" + p_samples[j] + ".temp", 'w')
            file_ext_temp.write('\t'.join(str(i) for i in res_fin))
            file_ext_temp.write('\n')
            file_ext_temp.close()

run("cat " + str(args['tempDir']) + "/*.temp > " + str(args['output']) + ".psi",shell=True)
run("rm -r " + str(args['tempDir']), shell = True)
