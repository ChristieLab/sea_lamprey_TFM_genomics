#!/bin/bash -l

#SBATCH --job-name=outlier_window_mRNA_blast_nr_rndII
#SBATCH -A beagle -p cpu -q standby
#SBATCH -t 4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --output="/scratch/bell/nbackens/Pmart_wgs/vcf_investigation/scripts/12_annot_all_metrics/blast_logs/outlier_window_blast_search_%j_%a.out"
#SBATCH --error="/scratch/bell/nbackens/Pmart_wgs/vcf_investigation/scripts/12_annot_all_metrics/blast_logs/outlier_window_blast_search_%j_%a.err"
#SBATCH --array=1-1854
echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
echo "SLURMTMPDIR="$SLURMTMPDIR
echo "SLURM_ARRAYID="$SLURM_ARRAYID
echo "SLURM_ARRAY_JOB_ID"=$SLURM_ARRAY_JOB_ID
echo "SLURM_ARRAY_TASK_ID"=$SLURM_ARRAY_TASK_ID
echo "working directory = "$SLURM_SUBMIT_DIR


module load biocontainers/default
module load blast


HOME="/scratch/bell/nbackens/Pmart_wgs/vcf_investigation/annot_all_metrics/" #
INPUT_DIR="/scratch/bell/nbackens/Pmart_wgs/vcf_investigation/annot_all_metrics/outlier_window_proteins_split/" #Input fasta directory for blast
OUTPUT_DIR="/scratch/bell/nbackens/Pmart_wgs/vcf_investigation/annot_all_metrics/outlier_window_mRNA_blastp_results/"
EVAL="1e-5" #e-value blast

mkdir -p ${OUTPUT_DIR} 

cd ${HOME}

RUN_ID=$(( $SLURM_ARRAY_TASK_ID + 1 ))
 
QUERY_FILE=$( ls ${INPUT_DIR} | sed -n ${RUN_ID}p )
QUERY_NAME="${QUERY_FILE%.*}"

echo "query name" $QUERY_NAME

QUERY="${INPUT_DIR}/${QUERY_FILE}"
OUTPUT="${OUTPUT_DIR}/${QUERY_NAME}.out"

blastp -query ${QUERY} -db nr -out ${OUTPUT} -evalue ${EVAL} -outfmt '6 std sseqid stitle' -max_hsps 50 -num_threads 4
 
date