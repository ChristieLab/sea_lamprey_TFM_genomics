#!/bin/bash -l

#SBATCH --job-name=Pmart_rFv3-2_whole
#SBATCH -A beagle
#SBATCH -t 336:00:00
#SBATCH -n 170
#SBATCH --output="./randomForest_chrom_v3-2_whole_100k.out"
#SBATCH --error="./randomForest_chrom_v3-2_whole_100k.err"


#Build reference index files
module load biocontainers/default
module load r
ulimit -s unlimited

Rscript rfv3.2-whole_VCF.r
