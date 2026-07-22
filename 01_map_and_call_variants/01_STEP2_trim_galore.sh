#!/bin/bash -l

#SBATCH --job-name=Pmart_trim
#SBATCH -A beagle
#SBATCH -t 220:00:00
#SBATCH -n 40
#SBATCH --output="/scratch/bell/nbackens/Pmart_wgs/logs/adapter_clip_logs/slurm-%A_%a.out"
#SBATCH --error="/scratch/bell/nbackens/Pmart_wgs/logs/adapter_clip_logs/slurm-%A_%a.err"
echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
echo "SLURMTMPDIR="$SLURMTMPDIR
echo "SLURM_ARRAYID="$SLURM_ARRAYID
echo "SLURM_ARRAY_JOB_ID"=$SLURM_ARRAY_JOB_ID
echo "SLURM_ARRAY_TASK_ID"=$SLURM_ARRAY_TASK_ID
echo "working directory = "$SLURM_SUBMIT_DIR


#REQUIRED INPUTS
MASTER_DIR=$1 # master directory containing data and where results will be written
#MASTER_FILE="sampleinfo.txt" # tab delimited master samples file, one row per sample, no headers
RAW_DATA=$2 #location of raw data directory
SUFFIX1=$3
SUFFIX2=$4

#LOAD MODULES
module load biocontainers/default
module load trim-galore/0.6.10
ulimit -s unlimited

#search all the fastq files from the "data" directory and generate the array
        index=$(( $SLURM_ARRAY_TASK_ID + 1 ))
        read1=$(ls ${RAW_DATA}/*${SUFFIX1} | sed -n ${index}p)
        read2=$(ls ${RAW_DATA}/*${SUFFIX2} | sed -n ${index}p)
        prefix=${read1%"$SUFFIX1"} # get file name prefix
        ID=${prefix#"${RAW_DATA}/"}

#STEP 2 - Trim galore
#navigate home
cd ${MASTER_DIR}

#generate a results directory
mkdir -p ./trimmed_reads

trim_galore -j 40 -o ${MASTER_DIR}/trimmed_reads --illumina --paired --retain_unpaired ${RAW_DATA}/${ID}${SUFFIX1} ${RAW_DATA}/${ID}${SUFFIX2}


