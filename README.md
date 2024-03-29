# Genomics tools with python
-
#### Author: Pierre Lesturgie
#### Example citation: Lesturgie, P. (2024). *Py-PopGen*. Available at: https://github.com/plesturgie/Py-PopGen. Last accessed XXX.
-

## Analyse and filter genomics data using commandline python scripts
 
 This repository provides python scripts to analyse and filter genomic data from VCF files. 
 
 For each script, all arguments are indicated using the --help argument.  
 
 Scripts run with python 3 (lastly tested with 3.10)

 To ensure that all script works, necessitates modules: 
 - pysam
 - pandas
 - tqdm
 - argparse
 - subprocess
 - itertools

### (1) Filtering SNPs with excess heterozygotes
Example with a threshold at 80% (output is  a vcf file named output.vcf)

	python HET.py --input input.vcf.gz --het 0.8 --output output 

List of arguments: 
- --input : input vcf file **must be compressed with bgzip and indexed with tabix (or uncompressed)**
- --het : threshold for heterozygote frequency (default is 0.8)
- --output : output prefix
- --chr : chromosome
- --write_discarded : write discarded SNPs (default is False)

### (2) Binning the vcf (for linkage disequilibrium)
There are two ways of binning: 
#### 1 - By regions, i.e., only keeps regions of <region> bp distant of <bin> from each other
Example with with regions of 100 bp distant of 100,000 bp to each other

	python binning.py --input input.vcf.gz --output binned_100k_regions_100bp.vcf.gz --bin 100000 --region 100

#### 2 - By SNP, i.e., only keeps SNPs distant of <bin> from each other
Example with with SNPs distant of 100,000 bp to each other

	python binning.py --input input.vcf.gz --output binned_100k.vcf.gz --bin 100000

List of arguments: 
- --input : input vcf file **must be compressed with bgzip and indexed with tabix (or uncompressed)**
- --bin : bin size
- --output : output name for vcf
- --region : region size

### (3) Polarize alleles of a vcf for a given outgroup individual
Example for a polarization of alleles for outgroup named "GN18033" (outputs: polarized.vcf.gz and polarized.der)

	python Polarize.py --input input_with_outgroup.vcf.gz --outgroup GN18033 --output polarized --tempDir tempola --write_derived_site

List of arguments: 
- --input : input vcf file 
- --outgroup : Name of the outgroup individual
- --output : output prefix
- --write_derived_site : Write a table with the number of derived alleles per individual per site
- --write_derived_ind : Write a table with the total number of derived alleles per individual **This option takes 5x more time to run**
- --tempDir : temporary directory
- --includeMiss : do not filter out missing data for vcf and .der file (default filters missing data)
- --include_outgroup : Include outgroup individual in the written VCF 

### (4) Site Frequency Spectrum 
There are two ways of computing the SFS: 
#### 1 - Accross all sites 
Remove --folded for unfolded sfs

	python SFS.py --input input.vcf.gz --folded --output folded.sfs

#### 2 - In sliding windows
Remove --folded for unfolded sfs
Example for contig SUPER_1 (--chr SUPER_1) in windows of 1,000,000 bp with a jump of 500,000 bp

	python SFS.py --input input.vcf.gz --folded --output sliding_folded.sfs --chr SUPER_1 --sliding_window --window 1000000 --jump 500000

List of arguments: 
- --input : input vcf file **must be compressed with bgzip and indexed with tabix (or uncompressed)**
- --folded : compute folded sfs (default is unfolded)
- --sliding_window : Computing in sliding window ? (Default is false)
- --jump : size of the jump (only in --sliding_window mode) 
- --window : size of the window (only in --sliding_window mode) 
- --chr CHR : chromosome
- --start : start of region
- --end : end of region
- --het : threshold for heterozygote frequency filtering (default is 1.0, i.e., no filtering)
- --output : output name
