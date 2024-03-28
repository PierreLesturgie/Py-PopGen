#!/bin/bash
# ********************************************************** #
# ---------------------------------------------------------- #
# This is a script to run in command line the python scripts #
# ---------------------------------------------------------- #
# ********************************************************** #

__author__      = "Pierre Lesturgie"
__email__       = "pierrelesturgie@outlook.fr"
__date__        = "2024-03-28"

### >> All scripts have a help: use --help in any to see additional arguments
### >> Tested with python 3.8

# --> 1. Filtering SNPs with excess heterozygotes
	# >> example with a threshold at 80% (output is  a vcf file named output.vcf)
echo "filtering sites with more than 80% of genotypes heterozygotes"
python HET.py --input input.vcf.gz -H 0.8 -O output 

# --> 2. Binning the vcf (for linkage disequilibrium)
	# >> There are two ways of binning: 
		# (1) by regions, i.e., only keeps regions of <region> bp distant of <bin> from each other
echo "binning regions of 100bp by 100k bins"
python binning.py --input input.vcf.gz --output binned_100k_regions_100bp.vcf.gz --bin 100000 --region 100

		# (2) by SNP, i.e., only keeps SNPs distant of <bin> from each other
echo "binning by 100k bins"
python binning.py --input input.vcf.gz --output binned_100k.vcf.gz --bin 100000

# --> 3. Polarize alleles of a vcf for a given outgroup individual
	# >> Example for a polarization of alleles for outgroup named "GN18033"
echo "Polarize alleles accordingly to given outgroup individual and writes a bgzipped vcf and a table with number of derived alleles per sample"
python Polarize.py --input input_with_outgroup.vcf.gz --outgroup GN18033 --output polarized_all_S1 -t tempola -D

# --> 4. Site Frequency Spectrum 
	# >> There are two ways of computing the SFS: 
			# (1) accross all sites (remove --folded for unfolded sfs)
echo "computing folded sfs including monomorphic sites in the output"
python SFS.py --input input.vcf.gz --folded --output folded.sfs

			# (2) in sliding windows (remove --folded for unfolded sfs)
echo "computing folded sfs in sliding windows (1Mb windows, 500kb jump) on the first contig and including monomorphic sites in the output"
python SFS.py --input input.vcf.gz --folded --output sliding_folded_S1.sfs --chr SUPER_1 --sliding_window --window 1000000 --jump 500000

