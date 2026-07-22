#!/bin/bash -l

#SBATCH --job-name=Pmart_map_sort_filt
#SBATCH -A beagle
#SBATCH -t 220:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --output="/scratch/bell/nbackens/Pmart_wgs/logs/logs_map_sort_filt/slurm-%A_%a.out"
#SBATCH --error="/scratch/bell/nbackens/Pmart_wgs/logs/logs_map_sort_filt/slurm-%A_%a.err"
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
MAPQUAL=$5
FILTERONLY=$6
STATSONLY=$7
RAM=$8
CPUs=${9}
RERUNSAMPLES=${10}
RERUNSAMPLEFILE=${11}

BOWTIE2_INDEXES=${REFPATH}

cd ${MASTER_DIR}

module load biocontainers/default
module load bwa/0.7.17
module load samtools/1.17
module load picard/2.26.10
module load bamutil/1.0.15
module load anaconda
module load use.own
module load conda-env/datamash
conda activate datamash

# only use the adapter clipped and paired reads from step 3
SUFFIX1="_1_val_1.fq.gz" #Read suffix forward
SUFFIX2="_2_val_2.fq.gz" #Read suffix reverse

if [[ ${RERUNSAMPLES} == "TRUE" ]]
then
	index=$(( $SLURM_ARRAY_TASK_ID + 1 ))
	ID=$(sed -n ${index}p ${MASTER_DIR}/${RERUNSAMPLEFILE})
	read1="${READ_DIR}/${ID}${SUFFIX1}"
	read2="${READ_DIR}/${ID}${SUFFIX2}"
else
	#search all the fastq files from the "data" directory and generate the array
	index=$(( $SLURM_ARRAY_TASK_ID + 1 ))
	read1=$(ls ${READ_DIR}/*${SUFFIX1} | sed -n ${index}p)
	read2=$(ls ${READ_DIR}/*${SUFFIX2} | sed -n ${index}p)
	prefix=${read1%"$SUFFIX1"} # get file name prefix
	ID=${prefix#"${READ_DIR}/"}
fi

echo "sample: ${ID}"
echo "${read1}"
echo "${read2}"

if [[ ${STATSONLY} == "TRUE" ]]
then
	echo "only calculating depth stats, assumes BAM files have already been generated"
	temp=$(samtools depth -aa ${MASTER_DIR}/bams_raw/${ID}_sorted.bam | cut -f 3 | datamash mean 1 median 1 min 1 max 1 sstdev 1)
	echo -e "${ID}\t${temp}" >> ${MASTER_DIR}/stats/bams_raw.stats.txt

	temp=$(samtools depth -aa ${MASTER_DIR}/bams_filtered/${ID}_sorted_minq${MAPQUAL}.bam | cut -f 3 | datamash mean 1 median 1 min 1 max 1 sstdev 1)
	echo -e "${ID}\t${temp}" >> ${MASTER_DIR}/stats/bams_minq${MAPQUAL}.stats.txt

	temp=$(samtools depth -aa ${MASTER_DIR}/bams_filtered/${ID}_sorted_minq${MAPQUAL}_dedup.bam | cut -f 3 | datamash mean 1 median 1 min 1 max 1 sstdev 1)
	echo -e "${ID}\t${temp}" >> ${MASTER_DIR}/stats/bams_minq${MAPQUAL}_dedup.stats.txt

	temp=$(samtools depth -aa ${MASTER_DIR}/bams_filtered/${ID}_sorted_minq${MAPQUAL}_dedup_overlapclipped.bam | cut -f 3 | datamash mean 1 median 1 min 1 max 1 sstdev 1)
	echo -e "${ID}\t${temp}" >> ${MASTER_DIR}/stats/bams_minq${MAPQUAL}_dedup_overlapclipped.stats.txt

	echo "depth calculation COMPLETE"


elif [[ ${RERUNSAMPLES} == "TRUE" ]]
then
	echo "starting initial read mapping and BAM sorting for samples listed in ${RERUNSAMPLEFILE}"
	# Map the paired-end reads
	# We ignore the reads that get orphaned during adapter clipping because that is typically a very small proportion of reads. If a large proportion of reads get orphaned (lose their mate so they become single-end), these can be mapped in a separate step and the resulting bam files merged with the paired-end mapped reads.
	rm -rf "${MASTER_DIR}/bams_raw/${ID}_sorted.bam"
	rm -rf "${MASTER_DIR}/bams_raw_${ID}_sorted.bam.temp*"
	bwa mem ${REFPATH}/${REFNAME} ${read1} ${read2} -t ${CPUs} | samtools view -buS | samtools sort -o ${MASTER_DIR}/bams_raw/${ID}_sorted.bam
	# get read depth stats
	temp=$(samtools depth -aa ${MASTER_DIR}/bams_raw/${ID}_sorted.bam | cut -f 3 | datamash mean 1 median 1 min 1 max 1 sstdev 1)
	echo -e "${ID}\t${temp}" >> ${MASTER_DIR}/stats/bams_raw.stats.txt

	echo "initial read mapping complete"

	echo "starting filtering"
	## Filter the mapped reads (to only retain reads with high mapping quality)
	# Filter bam files to remove poorly mapped reads (non-unique mappings and mappings with a quality score < 20)
	rm -rf "${MASTER_DIR}/bams_filtered/${ID}_sorted_minq${MAPQUAL}.bam"
	samtools view -h -q ${MAPQUAL} ${MASTER_DIR}/bams_raw/${ID}_sorted.bam | samtools view -buS -o ${MASTER_DIR}/bams_filtered/${ID}_sorted_minq${MAPQUAL}.bam
	# get read depth stats
	temp=$(samtools depth -aa ${MASTER_DIR}/bams_filtered/${ID}_sorted_minq${MAPQUAL}.bam | cut -f 3 | datamash mean 1 median 1 min 1 max 1 sstdev 1)
	echo -e "${ID}\t${temp}" >> ${MASTER_DIR}/stats/bams_minq${MAPQUAL}.stats.txt
	echo "finished filtering for mapping quality"

	echo "starting deduplicating"
	## Remove duplicates and print dupstat file
	rm -rf "${MASTER_DIR}/bams_filtered/${ID}_sorted_minq${MAPQUAL}_dedup.bam"
	picard MarkDuplicates I=./bams_filtered/${ID}_sorted_minq${MAPQUAL}.bam O=./bams_filtered/${ID}_sorted_minq${MAPQUAL}_dedup.bam M=./bams_filtered/${ID}_sorted_minq${MAPQUAL}_dupstat.txt VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true
	# get read depth stats
	temp=$(samtools depth -aa ${MASTER_DIR}/bams_filtered/${ID}_sorted_minq${MAPQUAL}_dedup.bam | cut -f 3 | datamash mean 1 median 1 min 1 max 1 sstdev 1)
	echo -e "${ID}\t${temp}" >> ${MASTER_DIR}/stats/bams_minq${MAPQUAL}_dedup.stats.txt
	echo "finished deduplicating"

	echo "starting to clip overlapping paired end reads"
	## Clip overlapping paired end reads (only necessary for paired-end data, so if you're only running se samples, you can comment this step out)
	rm -rf "${MASTER_DIR}/bams_filtered/${ID}_sorted_minq${MAPQUAL}_dedup_overlapclipped.bam"
	bam clipOverlap --in ./bams_filtered/${ID}_sorted_minq${MAPQUAL}_dedup.bam --out ./bams_filtered/${ID}_sorted_minq${MAPQUAL}_dedup_overlapclipped.bam --stats
	# get read depth stats
	temp=$(samtools depth -aa ${MASTER_DIR}/bams_filtered/${ID}_sorted_minq${MAPQUAL}_dedup_overlapclipped.bam | cut -f 3 | datamash mean 1 median 1 min 1 max 1 sstdev 1)
	echo -e "${ID}\t${temp}" >> ${MASTER_DIR}/stats/bams_minq${MAPQUAL}_dedup_overlapclipped.stats.txt
	echo "finished clipping overlapping paired end reads"
	echo "mapping and filtering COMPLETE"


elif [[ ${FILTERONLY} == "FALSE" ]]
then
	echo "starting initial read mapping and BAM sorting"
	# Map the paired-end reads
	# We ignore the reads that get orphaned during adapter clipping because that is typically a very small proportion of reads. If a large proportion of reads get orphaned (lose their mate so they become single-end), these can be mapped in a separate step and the resulting bam files merged with the paired-end mapped reads.
	if [ -f "${MASTER_DIR}/bams_raw/${ID}_sorted.bam" ]
	then
		echo "raw BAM file already exists, skipping"
	else
		bwa mem ${REFPATH}/${REFNAME} ${read1} ${read2} -t ${CPUs} | samtools view -buS | samtools sort -o ${MASTER_DIR}/bams_raw/${ID}_sorted.bam
		# get read depth stats
		temp=$(samtools depth -aa ${MASTER_DIR}/bams_raw/${ID}_sorted.bam | cut -f 3 | datamash mean 1 median 1 min 1 max 1 sstdev 1)
		echo -e "${ID}\t${temp}" >> ${MASTER_DIR}/stats/bams_raw.stats.txt
	fi
	echo "initial read mapping complete"

	echo "starting filtering"
	## Filter the mapped reads (to only retain reads with high mapping quality)
	# Filter bam files to remove poorly mapped reads (non-unique mappings and mappings with a quality score < 20)
	if [ -f "${MASTER_DIR}/bams_filtered/${ID}_sorted_minq${MAPQUAL}.bam" ]
	then
		echo "minq BAM file already exists, skipping"
	else
		samtools view -h -q ${MAPQUAL} ${MASTER_DIR}/bams_raw/${ID}_sorted.bam | samtools view -buS -o ${MASTER_DIR}/bams_filtered/${ID}_sorted_minq${MAPQUAL}.bam
		# get read depth stats
		temp=$(samtools depth -aa ${MASTER_DIR}/bams_filtered/${ID}_sorted_minq${MAPQUAL}.bam | cut -f 3 | datamash mean 1 median 1 min 1 max 1 sstdev 1)
		echo -e "${ID}\t${temp}" >> ${MASTER_DIR}/stats/bams_minq${MAPQUAL}.stats.txt
	fi
	echo "finished filtering for mapping quality"

	echo "starting deduplicating"
	## Remove duplicates and print dupstat file
	if [ -f "${MASTER_DIR}/bams_filtered/${ID}_sorted_minq${MAPQUAL}_dedup.bam" ]
	then
		echo "deduplicated BAM file already exists, skipping"
	else
		picard MarkDuplicates I=./bams_filtered/${ID}_sorted_minq${MAPQUAL}.bam O=./bams_filtered/${ID}_sorted_minq${MAPQUAL}_dedup.bam M=./bams_filtered/${ID}_sorted_minq${MAPQUAL}_dupstat.txt VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true
		# get read depth stats
		temp=$(samtools depth -aa ${MASTER_DIR}/bams_filtered/${ID}_sorted_minq${MAPQUAL}_dedup.bam | cut -f 3 | datamash mean 1 median 1 min 1 max 1 sstdev 1)
		echo -e "${ID}\t${temp}" >> ${MASTER_DIR}/stats/bams_minq${MAPQUAL}_dedup.stats.txt
	fi
	echo "finished deduplicating"

	echo "starting to clip overlapping paired end reads"
	## Clip overlapping paired end reads (only necessary for paired-end data, so if you're only running se samples, you can comment this step out)
	if [ -f "${MASTER_DIR}/bams_filtered/${ID}_sorted_minq${MAPQUAL}_dedup_overlapclipped.bam" ]
	then
		echo "overlap clipped BAM file already exists, skipping"
	else
		bam clipOverlap --in ./bams_filtered/${ID}_sorted_minq${MAPQUAL}_dedup.bam --out ./bams_filtered/${ID}_sorted_minq${MAPQUAL}_dedup_overlapclipped.bam --stats
		# get read depth stats
		temp=$(samtools depth -aa ${MASTER_DIR}/bams_filtered/${ID}_sorted_minq${MAPQUAL}_dedup_overlapclipped.bam | cut -f 3 | datamash mean 1 median 1 min 1 max 1 sstdev 1)
		echo -e "${ID}\t${temp}" >> ${MASTER_DIR}/stats/bams_minq${MAPQUAL}_dedup_overlapclipped.stats.txt
	fi
	echo "finished clipping overlapping paired end reads"
	echo "mapping and filtering COMPLETE"

else
	echo "skipping initial read mapping and BAM sorting"
	echo "starting filtering"
	## Filter the mapped reads (to only retain reads with high mapping quality)
	# Filter bam files to remove poorly mapped reads (non-unique mappings and mappings with a quality score < 20)
	if [ -f "${MASTER_DIR}/bams_filtered/${ID}_sorted_minq${MAPQUAL}.bam" ]
	then
		echo "minq BAM file already exists, skipping"
	else
		samtools view -h -q ${MAPQUAL} ${MASTER_DIR}/bams_raw/${ID}_sorted.bam | samtools view -buS -o ${MASTER_DIR}/bams_filtered/${ID}_sorted_minq${MAPQUAL}.bam
		# get read depth stats
		temp=$(samtools depth -aa ${MASTER_DIR}/bams_filtered/${ID}_sorted_minq${MAPQUAL}.bam | cut -f 3 | datamash mean 1 median 1 min 1 max 1 sstdev 1)
		echo -e "${ID}\t${temp}" >> ${MASTER_DIR}/stats/bams_minq${MAPQUAL}.stats.txt
	fi
	echo "finished filtering for mapping quality"

	echo "starting deduplicating"
	## Remove duplicates and print dupstat file
	if [ -f "${MASTER_DIR}/bams_filtered/${ID}_sorted_minq${MAPQUAL}_dedup.bam" ]
	then
		echo "deduplicated BAM file already exists, skipping"
	else
		java -jar /util/common/bioinformatics/picard/2.7.1/picard.jar MarkDuplicates I=./bams_filtered/${ID}_sorted_minq${MAPQUAL}.bam O=./bams_filtered/${ID}_sorted_minq${MAPQUAL}_dedup.bam M=./bams_filtered/${ID}_sorted_minq${MAPQUAL}_dupstat.txt VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true
		# get read depth stats
		temp=$(samtools depth -aa ${MASTER_DIR}/bams_filtered/${ID}_sorted_minq${MAPQUAL}_dedup.bam | cut -f 3 | datamash mean 1 median 1 min 1 max 1 sstdev 1)
		echo -e "${ID}\t${temp}" >> ${MASTER_DIR}/stats/bams_minq${MAPQUAL}_dedup.stats.txt
	fi
	echo "finished deduplicating"

	echo "starting to clip overlapping paired end reads"
	## Clip overlapping paired end reads (only necessary for paired-end data, so if you're only running se samples, you can comment this step out)
	if [ -f "${MASTER_DIR}/bams_filtered/${ID}_sorted_minq${MAPQUAL}_dedup_overlapclipped.bam" ]
	then
		echo "overlap clipped BAM file already exists, skipping"
	else
		bam clipOverlap --in ./bams_filtered/${ID}_sorted_minq${MAPQUAL}_dedup.bam --out ./bams_filtered/${ID}_sorted_minq${MAPQUAL}_dedup_overlapclipped.bam --stats
		# get read depth stats
		temp=$(samtools depth -aa ${MASTER_DIR}/bams_filtered/${ID}_sorted_minq${MAPQUAL}_dedup_overlapclipped.bam | cut -f 3 | datamash mean 1 median 1 min 1 max 1 sstdev 1)
		echo -e "${ID}\t${temp}" >> ${MASTER_DIR}/stats/bams_minq${MAPQUAL}_dedup_overlapclipped.stats.txt
	fi
	echo "finished clipping overlapping paired end reads"
	echo "filtering COMPLETE"
fi
