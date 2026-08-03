#!/bin/bash -l

#SBATCH --job-name=Pmart_WC_Fst
#SBATCH -A beagle
#SBATCH -t 20:00:00
#SBATCH -n 1
#SBATCH --output="./vcftools_Fst_calc.out"
#SBATCH --error="./vcftools_Fst_calc.err"


#Build reference index files
module load biocontainers/default
module load samtools/1.17
module load bcftools/1.17
module load vcftools/0.1.16
module load htslib/1.17
module load plink2
ulimit -s unlimited

HOME="/scratch/bell/nbackens/Pmart_wgs/vcf_investigation/" #Home directory
VCF="Pmart_filtered_CHROM_MLOCI0.8_MAF0.05_HWE0.05_QUAL30_MINDP10_MAXDP50_het_filt.recode.vcf.gz" #VCF file name
OUT="Pmart_filtered_CHROM_MLOCI0.8_MAF0.05_HWE0.05_QUAL30_MINDP10_MAXDP50_het_filt" #output file name prefix

cd ${HOME}

vcftools --gzvcf ${VCF} --weir-fst-pop ${HOME}/E.txt --weir-fst-pop ${HOME}/L.txt --out ./Fst/${OUT}_stdFst_EvL

vcftools --gzvcf ${VCF} --fst-window-size 50000 --fst-window-size 50000 --weir-fst-pop ${HOME}/E.txt --weir-fst-pop ${HOME}/L.txt --out ./Fst/${OUT}_50kbFst_EvL

vcftools --gzvcf ${VCF} --fst-window-size 100000 --fst-window-size 100000 --weir-fst-pop ${HOME}/E.txt --weir-fst-pop ${HOME}/L.txt --out ./Fst/${OUT}_100kbFst_EvL
