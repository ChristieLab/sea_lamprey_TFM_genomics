#!/bin/bash -l

#SBATCH --job-name=Pmart_index_ref
#SBATCH -A beagle
#SBATCH -t 220:00:00
#SBATCH -n 40
#SBATCH --output="/scratch/bell/nbackens/Pmart_wgs/logs/index_logs/Pmart_index.out"
#SBATCH --error="/scratch/bell/nbackens/Pmart_wgs/logs/index_logs/Pmart_index.err"

echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
echo "SLURMTMPDIR="$SLURMTMPDIR
echo "SLURM_ARRAYID="$SLURM_ARRAYID
echo "SLURM_ARRAY_JOB_ID"=$SLURM_ARRAY_JOB_ID
echo "SLURM_ARRAY_TASK_ID"=$SLURM_ARRAY_TASK_ID
echo "working directory = "$SLURM_SUBMIT_DIR


#Build reference index files
module load biocontainers/default
module load samtools/1.17
module load picard/2.26.10
module load bowtie2/2.5.1
module load bwa/0.7.17
ulimit -s unlimited


MASTER_DIR=$1 
REF=$2   # This is a fasta file with the reference genome sequence we will map to 
REFBASENAME="${REF%.fna*}"

cd ${MASTER_DIR}

samtools faidx ${REF}

picard CreateSequenceDictionary R=${REF} O=${REFBASENAME}'.dict'

bowtie2-build ${REF} ${REFBASENAME}

bwa index ${REF}