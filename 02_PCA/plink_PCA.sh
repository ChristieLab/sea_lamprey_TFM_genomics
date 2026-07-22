#!/bin/bash -l

#SBATCH --job-name=Pmart_plink_PCA
#SBATCH -A beagle
#SBATCH -t 20:00:00
#SBATCH -n 40
#SBATCH --output="./plinkPCAII.out"
#SBATCH --error="./plinkPCAII.err"


#Build reference index files
module load biocontainers/default
module load samtools/1.17
module load bcftools/1.17
module load vcftools/0.1.16
module load htslib/1.17
module load plink2
ulimit -s unlimited

HOME="/scratch/bell/nbackens/Pmart_wgs/tempVI/" #Home directory
VCF="Pmart_filtered_MLOCI0.8_MAF0.05_HWE0.05_QUAL30_MINDP10_MAXDP50.vcf.gz" #VCF file name
OUT="Pmart_filtered_MLOCI0.8_MAF0.05_HWE0.05_QUAL30_MINDP10_MAXDP50" #output file name prefix

cd ${HOME}

#plink2 --vcf ${VCF} --make-bed --allow-extra-chr --out ${OUT}

plink2 --bfile ${OUT} --pca 223 --allow-extra-chr --out ${OUT}_plinkPCA