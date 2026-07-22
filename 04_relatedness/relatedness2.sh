#!/bin/bash -l

#SBATCH --job-name=Pmar_relatedness
#SBATCH -A beagle
#SBATCH -t 180:00:00
#SBATCH -n 1
#SBATCH --output="./Pmart_relatedness2_het_filt.out"
#SBATCH --error="./Pmart_relatedness2_het_filt.err"


#Build reference index files
module load biocontainers/default
module load samtools/1.17
module load bcftools/1.17
module load vcftools/0.1.16
module load htslib/1.17
ulimit -s unlimited

#INPUTS
HOME="/scratch/bell/nbackens/Pmart_wgs/vcf_investigation/relatedness_het_filt" #Home directory
VCF="Pmart_filtered_CHROM_MLOCI0.8_MAF0.05_HWE0.05_QUAL30_MINDP10_MAXDP50_het_filt.recode.vcf.gz" #VCF file name
OUT="Pmart_filtered_CHROM_het_filt" #output prefix, for me it's the gene name, chromosome, and coordinates 
                                        #used in various locations throughout

#GO HOME 
cd ${HOME}

#STEP 1 - Calculate relatedness
vcftools --gzvcf ${VCF} --relatedness2 --out ${OUT}
