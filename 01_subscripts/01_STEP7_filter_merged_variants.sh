#!/bin/bash -l

#SBATCH --job-name=Pmart_bcftools_filter
#SBATCH -A beagle
#SBATCH -t 220:00:00
#SBATCH -n 40
#SBATCH --output="/scratch/bell/nbackens/Pmart_wgs/logs/vcftools_filter_logsII/vcftools_filter_variants.out"
#SBATCH --error="/scratch/bell/nbackens/Pmart_wgs/logs/vcftools_filter_logsII/vcftools_filter_variants.err"

echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
echo "SLURMTMPDIR="$SLURMTMPDIR
echo "SLURM_ARRAYID="$SLURM_ARRAYID
echo "SLURM_ARRAY_JOB_ID"=$SLURM_ARRAY_JOB_ID
echo "SLURM_ARRAY_TASK_ID"=$SLURM_ARRAY_TASK_ID
echo "working directory = "$SLURM_SUBMIT_DIR

MASTER_DIR=$1
CPUs=$2
CLEANUP=$3
ID=$4
MIND=$5
MAF=$6
QUAL=$7
MIN_DEPTH=$8
MAX_DEPTH=$9
MAPQUAL=$10


#Build reference index files
module load biocontainers/default
module load samtools/1.17
module load bcftools/1.17
module load vcftools/0.1.16

ulimit -s unlimited

cd ${MASTER_DIR}
cd ${MASTER_DIR}/combined_variantsIII

#Filter for quality and depth

bcftools view --threads ${CPUs} -m2 -M2 -v snps -Oz -o /scratch/bell/nbackens/Pmart_wgs/combined_variantsIII/Pmart_bamMinQ20_ALL_VARIANTS_reheader_biallelicSNPs.vcf.gz /scratch/bell/nbackens/Pmart_wgs/combined_variantsIII/Pmart_bamMinQ20_ALL_VARIANTS_reheader.vcf.gz

vcftools --gzvcf /scratch/bell/nbackens/Pmart_wgs/combined_variantsIII/Pmart_bamMinQ20_ALL_VARIANTS_reheader_biallelicSNPs.vcf.gz --max-missing $MIND --recode --stdout | gzip -c > temp_${ID}_MIND${MIND}.vcf.gz

vcftools --gzvcf temp_${ID}_MIND${MIND}.vcf.gz --maf $MAF --recode --stdout | gzip -c > temp_${ID}_MIND${MIND}_MAF${MAF}.vcf.gz

vcftools --gzvcf temp_${ID}_MIND${MIND}_MAF${MAF}.vcf.gz --minQ $QUAL  --recode --stdout | gzip -c > temp_${ID}_MIND${MIND}_MAF${MAF}_QUAL${QUAL}.vcf.gz

vcftools --gzvcf temp_${ID}_MIND${MIND}_MAF${MAF}_QUAL${QUAL}.vcf.gz \
--minDP ${MIN_DEPTH} --maxDP ${MAX_DEPTH} \
--recode --stdout | gzip -c > ${MASTER_DIR}/filtered_biallelic_snpsII/${ID}_filtered_MIND${MIND}_MAF${MAF}_QUAL${QUAL}_MINDP${MIN_DEPTH}_MAXDP${MAX_DEPTH}.vcf.gz

if [ "$my_variable" = "TRUE" ]; then
    rm ./temp*
else
    mv ./temp* ${MASTER_DIR}/filtered_biallelic_snpsII/
fi
