#!/bin/bash -l

#SBATCH --job-name=Pmart_merge_variants
#SBATCH -A beagle
#SBATCH -t 220:00:00
#SBATCH -n 40
#SBATCH --output="/scratch/bell/nbackens/Pmart_wgs/logs/bcftools_merge_logsIII/bcftools_merge_variants.out"
#SBATCH --error="/scratch/bell/nbackens/Pmart_wgs/logs/bcftools_merge_logsIII/bcftools_merge_variants.err"

echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
echo "SLURMTMPDIR="$SLURMTMPDIR
echo "SLURM_ARRAYID="$SLURM_ARRAYID
echo "SLURM_ARRAY_JOB_ID"=$SLURM_ARRAY_JOB_ID
echo "SLURM_ARRAY_TASK_ID"=$SLURM_ARRAY_TASK_ID
echo "working directory = "$SLURM_SUBMIT_DIR

MASTER_DIR=$1
MAPQUAL=$2
CPUs=$3


#Build reference index files
module load biocontainers/default
module load samtools/1.17
module load bcftools/1.17

ulimit -s unlimited

cd ${MASTER_DIR}
cd ${MASTER_DIR}/variants_rawII

#create list of prefixes
SUFFIX="_bamMinQ${MAPQUAL}_multiallelic.vcf.gz"
ls *${SUFFIX} > temp
#sed -i "s/${SUFFIX}//g" "temp"
#FILELIST=`cat temp`
#rm temp


#merge variants
bcftools merge --threads ${CPUs} -l temp -Oz -o Pmart_bamMinQ${MAPQUAL}_ALL_VARIANTS.vcf.gz

#cmoveto new directory

mv Pmart_bamMinQ${MAPQUAL}_ALL_VARIANTS.vcf.gz ${MASTER_DIR}/combined_variantsIII/Pmart_bamMinQ${MAPQUAL}_biallelic_SNPs.vcf.gz


#for f in $FILELIST
#do
##filter for biallelic SNPs
#bcftools view --threads ${CPUs} -m2 -M2 -v snps -Oz ${f}${SUFFIX} > ${MASTER_DIR}/biallelic_snpsII/${f}_bamMinQ${MAPQUAL}_biallelic_SNPs_raw.vcf.gz
#index
#bcftools index --threads ${CPUs} ${MASTER_DIR}/biallelic_snpsII/${f}_bamMinQ${MAPQUAL}_biallelic_SNPs_raw.vcf.gz
#done






