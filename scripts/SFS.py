#!/usr/bin/python

__author__      = "Pierre Lesturgie"
__version__     = "4.1"
__email__       = "pierrelesturgie@outlook.fr"
__date__        = "2024-03-07"

from pysam import VariantFile
from pandas import DataFrame
from tqdm import tqdm
from argparse import ArgumentParser
from subprocess import run

# Construct the argument parser
ap = ArgumentParser(description='Compute Unfolded and/or Folded Site Frequency Spectrum from vcf')

# Add the arguments to the parser
ap.add_argument("-i", "--input", required=True,
                help="input vcf file must be compressed with bgzip and indexed with tabix (or uncompressed)")
ap.add_argument("-f", "--folded", required=False,
                help="compute Folded sfs (default is unfolded)",action='store_true', default=False)
ap.add_argument("-s", "--sliding_window", required=False,
                help="Computing in sliding window ? (Default is false)",action='store_true', default=False)
ap.add_argument("-j", "--jump", required=False,
                help="size of the jump",type=int,default=None)
ap.add_argument("-w", "--window", required=False,
                help="size of the window",type=int,default=None)
ap.add_argument("-c", "--chr", required=False,
                help="chromosome",type=str,default=None)
ap.add_argument("-b", "--start", required=False,
                help="start of region",type=int,default=None)
ap.add_argument("-e", "--end", required=False,
                help="end of region",type=int,default=None)
ap.add_argument("-H", "--het", required=False,
                help="threshold for heterozygote frequency (default is 1.0, i.e., no filtering)",type=float, default=1.0)
ap.add_argument("-O", "--output", required=False,
                help="output name",type=str, default="out.sfs")

args = vars(ap.parse_args())

#### The sliding window is now correct : before, I took the end of the vcf file in bash as the last position
#### but it was dumb for multi-contigs datasets
#### Now I just extract the length from the header in bash using temporary files
#### Computation time is overall smaller: the code has changed a lot (hence v4)
#### Sliding window mode now works without specifying a contig: Automatic sep of contigs. 
#### But will also loop on ghost contigs, so will output a lot of non important warnings which indicates 
#### no output for ghosts (otherwise, very embarrassing amount of output data). 

#### v.4.1: two modifications: (1) VCF-compatible with output from DEEPVARIANT; (2) filtering out missing data


if args['sliding_window'] == True: 
    if args['window'] is None or args['jump'] is None :
        sys.exit("Need a window and jump value in sliding_window mode")

if args['end'] is not None: 
    if args['chr'] is None or args['start'] is None:
        sys.exit("Need a chromosome, start and end value in subset region mode")
        
if args['start'] is not None: 
    if args['chr'] is None or args['end'] is None:
        sys.exit("Need a chromosome, start and end value in subset region mode")
    
    
############## Definition of functions used
#### Function 1: counts number of derived allele
def derived_allele_count(rec,het_threshold):
    alleles = []
    het = 0
    filter_out=False
    for s in rec.samples.values() :
        if s['GT'][0] is not None and s['GT'][1] is not None: 
            alleles.append(s['GT'][0])
            alleles.append(s['GT'][1])
            if s['GT']==(0,1) or s['GT']==(1,0):
                het = het + 1
        else: filter_out=True

    if filter_out == False:
        if len([item for item in alleles if str(item) != 'None']) == len(alleles):
            if (het / (len(alleles) / 2)) <= het_threshold:
                return(sum(alleles)) ## Now only rturns IF and only IF the site passes the het and missing data filter

def derived_allele_count_old(rec,het_threshold):
    alleles = []
    het = 0
    for s in rec.samples.values() :
        alleles.append(s['GT'][0])
        alleles.append(s['GT'][1])
        if s['GT']==(0,1) or s['GT']==(1,0):
            het = het + 1
    if len([item for item in alleles if str(item) != 'None']) == len(alleles):
        if (het / (len(alleles) / 2)) <= het_threshold:
            return(sum(alleles)) ## Now only rturns IF and only IF the site passes the het filter

#### Function 2: computes unfolded SFS
def unfolded_SFS(vcf_in,het_threshold=0.8,REGION=None,START=None,STOP=None, sliding_windows=False):
    
    n_samples = len(list((vcf_in.header.samples)))
    n_derived_alleles = []
    
    if sliding_windows == True:
        for rec in vcf_in.fetch(start=START, end=STOP, region=REGION):
            n_derived_alleles.append(derived_allele_count(rec,het_threshold))      
    else:
        for rec in tqdm(vcf_in.fetch(start=START, end=STOP, region=REGION)):
            n_derived_alleles.append(derived_allele_count(rec,het_threshold))  
    

    df_nda = DataFrame({"Number_of_Derived_alleles":n_derived_alleles})

    # computing the unfolded sfs
    derived_allele_SFS = df_nda["Number_of_Derived_alleles"].value_counts().sort_index()

    temp_unfolded=[]
    for i in range(n_samples*2 + 1):
        if i in derived_allele_SFS.index:
            temp_unfolded.append(derived_allele_SFS[i])
        else:
            temp_unfolded.append(0)
    return(temp_unfolded)

#### Function 2: Folds unfolded SFS
def fold_SFS(temp_unfolded,n_samples):
    foldedSFS = [temp_unfolded.pop(0) + temp_unfolded.pop(-1)]
                
    for k in range(n_samples-1):
        foldedSFS.append(temp_unfolded[k] + temp_unfolded[-k -1])
                
    foldedSFS.append(temp_unfolded[n_samples-1])
    return(foldedSFS)

def sliding_function(pos,window,slide,het,CONTIG,FOLD):
    print("<<< Computing sliding windows SFS for contig " + str(CONTIG) + "...", flush=True)
    for i in tqdm(range(1,pos,window)):
        if i !=1:
            for w in slide:
                temp_unfolded = unfolded_SFS(vcf_in,het_threshold=het,REGION=CONTIG,START=temp + w,STOP=i + w,sliding_windows=True)
                d = (temp + w + i + w)/2
                
                if FOLD == True:
                    foldedSFS = fold_SFS(temp_unfolded,n_samples)
                    res =  DataFrame(foldedSFS).transpose()
                    suffix = "_folded"
                else : 
                    res =  DataFrame(temp_unfolded).transpose()
                    suffix = "_unfolded"
                
                res.insert(0, "pos", round(d))
                res.insert(0, "chrom", str(CONTIG))
                
                res.to_csv("temp_dir_sfs/" + CONTIG + "_"+ str(round(d)) + suffix + ".sfs",index=False, sep=" ",header=False)
                
            temp = i
        else:
            temp = i
    print("...Done >>>", flush=True)
    


######## Beginning of the "real" code
# Reading input VCF file
vcf_in = VariantFile(str(args['input']))
n_samples = len(list((vcf_in.header.samples)))

if args['sliding_window'] == True:
    
    run(["mkdir", "temp_dir_sfs"])
    
    print("<<< Extracting chromosome length info...", flush=True)
    run("zgrep contig " + str(args['input']) + "> temp_dir_sfs/_temp",shell=True)
    file_temp_loc = open('temp_dir_sfs/_temp', 'r')
    Lines = file_temp_loc.readlines()
    
    CONTIGS, LENGTH = [], []
    for i in Lines: 
        temp = i.split(",")[0]
        CONTIGS.append(temp.split("=")[2])
        temp = i.split(",")[1]
        LENGTH.append(temp.split("=")[1])
    
    LENGTH = [
    item.replace('\r', '').replace('>\n', '').replace('.\n', '')
    for item in LENGTH
    ]
    
    run("rm temp_dir_sfs/_temp",shell=True)
    
    print("...Done >>>", flush=True)
    
    slide = range(0, args['window'], args['jump'])
    
    if args['chr'] is not None :  
        pos = int(LENGTH[CONTIGS.index(str(args['chr']))])
        sliding_function(pos,window=args['window'],slide=slide,het=args['het'],CONTIG=args['chr'],FOLD=args['folded'])
    
    else : 
        for POSITION in LENGTH:
            CONTIG_temp = CONTIGS[LENGTH.index(POSITION)]
            
            contig_occ=[]
            for rec in vcf_in.fetch(str(CONTIG_temp)):
                contig_occ.append(rec.contig)
                break
                
            if len(contig_occ) > 0:
                pos = int(POSITION)
                sliding_function(pos,window=args['window'],slide=slide,het=args['het'],
                                 CONTIG=CONTIG_temp,FOLD=args['folded'])
            else:
                print("WARNING >>> No data detected for contig " + CONTIG_temp + "... Not writting output", flush=True)
            
    run(["cat temp_dir_sfs/*.sfs > ./sliding_window_" + args['output']],shell=True)
    run(["rm -r temp_dir_sfs"],shell=True)


if args['sliding_window'] == False:
    temp_unfolded = unfolded_SFS(vcf_in,het_threshold=args['het'],REGION=args['chr'],START=args['start'],STOP=args['end'])
    if args['folded'] == False:
        derived_allele_SFS = DataFrame({"Number_of_Derived_alleles":temp_unfolded})
        DataFrame(derived_allele_SFS).transpose().to_csv(args['output'],index=False, sep=" ")

    if args['folded'] == True:
        foldedSFS = fold_SFS(temp_unfolded,n_samples)
        DataFrame(foldedSFS).transpose().to_csv(args['output'],index=False, sep=" ")

