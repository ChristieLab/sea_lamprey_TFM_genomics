#!/bin/bash -l

#SBATCH --job-name=Pmart_bcftools
#SBATCH -A beagle
#SBATCH -t 220:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --output="/scratch/bell/nbackens/Pmart_wgs/logs/bcftools_logs/slurm-%A_%a.out"
#SBATCH --error="/scratch/bell/nbackens/Pmart_wgs/logs/bcftools_logs/slurm-%A_%a.err"
echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
echo "SLURMTMPDIR="$SLURMTMPDIR
echo "SLURM_ARRAYID="$SLURM_ARRAYID
echo "SLURM_ARRAY_JOB_ID"=$SLURM_ARRAY_JOB_ID
echo "SLURM_ARRAY_TASK_ID"=$SLURM_ARRAY_TASK_ID
echo "working directory = "$SLURM_SUBMIT_DIR


MASTER_DIR=$1
READ_DIR=$2
REFPATH=$3
REFNAME=$4
REFBASENAME="${REFNAME%.fna*}"
BAMS_FILTERED=$5
MAPQUAL=$6
CPUs=${7}
C=${8}
MAX_DEPTH=${9}

cd ${MASTER_DIR}

module load biocontainers/default
module load samtools/1.17
module load bcftools/1.17

# only use the adapter clipped and paired reads from step 3 to generate array
SUFFIX="_sorted_minq${MAPQUAL}_dedup_overlapclipped.bam" #Read suffix for filtered bam file


	#search all the bam files from the alignment directory and generate the array
	index=$(( $SLURM_ARRAY_TASK_ID + 1 ))
  bam_file=$(ls ${BAMS_FILTERED}/*${SUFFIX} | sed -n ${index}p)
  prefix=${bam_file%"$SUFFIX"}
  ID=${prefix#"${BAMS_FILTERED}/"}

echo "The identifier for task ${SLURM_ARRAY_TASK_ID} is: ${ID}"

#go to directory of filtered bams
cd ${MASTER_DIR}
cd ${BAMS_FILTERED}

#index overlap clipped, deduplicated, minq${MAPQUAL} bam file
samtools index -@ ${CPUs} ${ID}_sorted_minq${MAPQUAL}_dedup_overlapclipped.bam

#summarize the alignment of reads to a reference sequence at each position |\
#call SNPs, reporting only variant sites and multiallelic variants
bcftools mpileup --threads ${CPUs} -C ${C} -d ${MAX_DEPTH} -f ${REFPATH}/${REFNAME} ${ID}_sorted_minq${MAPQUAL}_dedup_overlapclipped.bam -Ou | \
bcftools call --threads ${CPUs} -mv -Oz > ${MASTER_DIR}/variants_raw/${ID}_bamMinQ${MAPQUAL}_multiallelic.vcf.gz

#index output
cd ${MASTER_DIR}/variants_raw/
bcftools index ${ID}_bamMinQ${MAPQUAL}_multiallelic.vcf.gz




