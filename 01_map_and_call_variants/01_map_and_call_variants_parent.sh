#!/bin/bash
############################################################
# ABOUT                                                    #
############################################################
#Petromyzon marinus map and call variants pipeline
#Modified from Dan MacGuigan SNP calling pipeline (INSERT GITHUB)
#author: Nate Backenstose
#date: 06/12/2024


############################################################
# INPUTS                                                   #
############################################################
#REQUIRED INPUTS
MASTER_DIR="/scratch/bell/nbackens/Pmart_wgs" # master directory containing data and where results will be written
#MASTER_FILE="sampleinfo.txt" # tab delimited master samples file, one row per sample, no headers
RAW_DATA="/scratch/bell/nbackens/Pmart_wgs/Pmart_raw_reads" #location of raw data directories
SUFFIX1="_1.fq.gz" #Raw read suffix forward
SUFFIX2="_2.fq.gz" #Raw read suffix reverse
NSAMPLES=225 #number of samples
REFPATH="/scratch/bell/nbackens/Pmart_wgs/reference/" # full path to reference genome
REFNAME="GCF_010993605.1_kPetMar1.pri_genomic.fna" # name of reference genome file in REFPATH, must end in .fasta
READ_DIR="/scratch/bell/nbackens/Pmart_wgs/trimmed_reads/" #Location of trimmed reads (fill in after step 2)
ID="Pmart" #Name 

#INDEX OPTIONS
#N/A

#ALIGNMENT OPTIONS
RAM=200 # Gb of RAM for mapping and filtering
CPUs=20 # number of CPUs for mapping and filtering
#MAPPINGPRESET="very-sensitive" # bowtie2 mapping preset, we recommend "very-sensitive" http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#preset-options-in---end-to-end-mode
MAPQUAL=20 # mapping quality cutoff, reads with lower mapq score will be discarded
FILTERONLY="FALSE" # only sort and filter, TRUE or FALSE? Useful if you have already generated raw BAM files in step 6 and want to apply a new MAPQUAL filter
STATSONLY="FALSE" # only generate depth stats, do no map reads to reference, TRUE or FALSE. Assumes all BAM files already exist
RERUNSAMPLES="FALSE" # rerun subset of samples? TRUE or FALSE. This will delete any BAM files that already exist for these samples. Need to generate new stats files after rerunning samples (use STATSONLY="TRUE")
RERUNSAMPLESFILE="rerun.txt" # file of sample IDs to rerun, one line per sample, should be within MASTER_DIR

#SNP calling options (BCFtools)
BAMS_FILTERED="${MASTER_DIR}/bams_filtered" #directory containing filtered bam files [you shouldn't need to change this]
#CPUs=adjust alignment options setting above or insert new value here and comment out line above
C=50 #mpileup - Coefficient for downgrading mapping quality for reads containing excessive mismatches. 
        #Given a read with a phred-scaled probability q of being generated from the mapped position, the new mapping quality is about sqrt((INT-q)/INT)*INT. 
        #A zero value (the default) disables this functionality.
#MAX_DEPTH=250 #mpileup - At a position, read maximally INT reads per input file. [250 is default]

#FILTERING OPTIONS 
MIND="0.8" #missingness per sample
MAF="0.05" #minor allele freqency value
MIN_DEPTH="10"
MAX_DEPTH="50"
QUAL="30"
CLEANUP="FALSE" #Remove temporary vcf files during filtering


############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo ""
   echo "WGS map-sort-filter pipeline"
   echo "Here we can describe each step in the pipeline"
   echo "Step 1: fastqc"
   echo "Step 2: adapter clipping with trim_galore"
   echo "Step 3: index the reference fasta"
   echo "Step 4: map sort and filter reads and collect stats"
   echo "Step 5: call variants with bcftools"
   echo "Step 6: merge raw variants with bcftools"
   echo "Step 7: filter variants with vcftools"
   echo "Syntax: scriptTemplate [-s|h]"
   echo "options:"
   echo "s    Which step would you like to run? Enter number 1 through N"
   echo "h    Print this Help"
   echo
}

############################################################
# Process the input options. Add options as needed.        #
############################################################
# Get the options
while getopts ":hs:" option; do
   case $option in
      h) # display Help
         Help
         exit;;
      s) # run a step in the pipeline
         step=$OPTARG
        echo "running step ${step}";;
      \?) # Invalid option
         echo "Error: Invalid option"
        Help
        exit;;
   esac
done


############################################################
# STEPS                                                    #
############################################################
#STEP 1 - Fastqc
if [[ "${step}" == 1 ]]; then
        mkdir -p ${MASTER_DIR}/logs
        mkdir -p ${MASTER_DIR}/logs/fastqc_logs
        echo "submitting the following job:"
        echo "sbatch --array=0-$(( NSAMPLES * 2 - 1 )) ./01_subscripts/01_STEP1_fastqc.sh ${MASTER_DIR} ${RAW_DATA} ${SUFFIX1} ${SUFFIX2}"
        sbatch --array=0-$(( NSAMPLES * 2 - 1 )) ./01_subscripts/01_STEP1_fastqc.sh ${MASTER_DIR} ${RAW_DATA} ${SUFFIX1} ${SUFFIX2}
        echo "log files saved in ${MASTER_DIR}/logs/fastqc_logs"
fi


#STEP 2 - Adapter clipping
if [[ "${step}" == 2 ]]; then
        echo "submitting the following job:"
        echo "sbatch --array=0-$(( NSAMPLES - 1 )) ./01_subscripts/0_STEP2_trim_galore.sh ${MASTER_DIR} ${RAW_DATA} ${SUFFIX1} ${SUFFIX2}"
        mkdir -p ${MASTER_DIR}/logs/adapter_clip_logs
        sbatch --array=0-$(( NSAMPLES - 1 )) ./01_subscripts/01_STEP2_trim_galore.sh ${MASTER_DIR} ${RAW_DATA} ${SUFFIX1} ${SUFFIX2} 
        echo "log files saved in ${MASTER_DIR}/adapter_clip_logs"
fi


#STEP 3 - INDEX REFERENCE
if [[ "${step}" == 3 ]]; then
        echo "submitting the following job:"
        echo "sbatch ./01_subscripts/01_STEP3_index.sh ${MASTER_DIR} ${REFPATH}/${REFNAME}"
        mkdir -p ${MASTER_DIR}/logs/index_logs
        sbatch ./01_subscripts/01_STEP3_index.sh ${MASTER_DIR} ${REFPATH}/${REFNAME}
        echo "log files saved in ${MASTER_DIR}/logs/index_logs"
fi


#STEP 4 - Map reads to ref sort and filter
if [[ "${step}" == 4 ]]; then
	echo "submitting the following job:"
  if [[ ${RERUNSAMPLES} == "TRUE" ]]
  then
    NRERUNSAMPLES=$(wc -l ${MASTER_DIR}/${RERUNSAMPLESFILE} | cut -f1 -d" ")
  	echo "sbatch --array=0-$(( NRERUNSAMPLES - 1 )) ./01_subscripts/01_STEP4_map_sort_filter.sh ${MASTER_DIR} ${READ_DIR} ${REFPATH} ${REFNAME} ${MAPQUAL} ${FILTERONLY} ${STATSONLY} ${RAM} ${CPUs} ${RERUNSAMPLES} ${RERUNSAMPLESFILE}"
  else
    echo "sbatch --array=0-$(( NSAMPLES - 1 )) ./01_subscripts/01_STEP4_map_sort_filter.sh ${MASTER_DIR} ${READ_DIR} ${REFPATH} ${REFNAME} ${MAPQUAL} ${FILTERONLY} ${STATSONLY} ${RAM} ${CPUs} ${RERUNSAMPLES} ${RERUNSAMPLESFILE}"
  fi
  mkdir -p ${MASTER_DIR}/logs/logs_map_sort_filt
	mkdir -p ${MASTER_DIR}/bams_raw
	mkdir -p ${MASTER_DIR}/bams_filtered
  if [[ ${RERUNSAMPLES} == "TRUE" ]]
  then
    NRERUNSAMPLES=$(wc -l ${MASTER_DIR}/${RERUNSAMPLESFILE} | cut -f1 -d" ")
  	sbatch --array=0-$(( NRERUNSAMPLES - 1 )) ./01_subscripts/01_STEP4_map_sort_filter.sh ${MASTER_DIR} ${READ_DIR} ${REFPATH} ${REFNAME}  ${MAPQUAL} ${FILTERONLY} ${STATSONLY} ${RAM} ${CPUs} ${RERUNSAMPLES} ${RERUNSAMPLESFILE}
  else
    rm ${MASTER_DIR}/stats/bams_raw.stats.txt
    rm ${MASTER_DIR}/stats/bams_minq${MAPQUAL}.stats.txt
    rm ${MASTER_DIR}/stats/bams_minq${MAPQUAL}_dedup.stats.txt
    rm ${MASTER_DIR}/stats/bams_minq${MAPQUAL}_dedup_overlapclipped.stats.txt
    sbatch --array=0-$(( NSAMPLES - 1 )) ./01_subscripts/01_STEP4_map_sort_filter.sh ${MASTER_DIR} ${READ_DIR} ${REFPATH} ${REFNAME} ${MAPQUAL} ${FILTERONLY} ${STATSONLY} ${RAM} ${CPUs} ${RERUNSAMPLES} ${RERUNSAMPLESFILE}
  fi
	echo "raw mapping results in ./bams_raw"
	echo "filtered and sorted mapping results in ./bams_filtered"
fi

#STEP 5 - SNP calling with bcftools
if [[ "${step}" == 5 ]]; then
        echo "submitting the following job:"
        echo "sbatch --array=0-$(( NSAMPLES - 1 )) ./01_subscripts/01_STEP5_bcftools.sh ${MASTER_DIR} ${READ_DIR} ${REFPATH} ${REFNAME} ${BAMS_FILTERED} ${MAPQUAL} ${CPUs} ${C} ${MAX_DEPTH}"
        mkdir -p ${MASTER_DIR}/logs/bcftools_logs
        mkdir -p ${MASTER_DIR}/variants_raw
        sbatch --array=0-$(( NSAMPLES - 1 )) ./01_subscripts/01_STEP5_bcftools.sh ${MASTER_DIR} ${READ_DIR} ${REFPATH} ${REFNAME} ${BAMS_FILTERED} ${MAPQUAL}  ${CPUs} ${C} ${MAX_DEPTH}
        echo "log files saved in ${MASTER_DIR}/logs/bcftools_logs"
fi


#STEP 6 - merge raw variants
if [[ "${step}" == 6 ]]; then
        echo "submitting the following job:"
        echo "sbatch ./01_subscripts/01_STEP6_merge_variants2.sh ${MASTER_DIR} ${MAPQUAL} ${CPUs}"
        mkdir -p ${MASTER_DIR}/logs/bcftools_merge_logsIII
        #mkdir -p ${MASTER_DIR}/biallelic_snpsIII
        mkdir -p ${MASTER_DIR}/combined_variantsIII
        sbatch ./01_subscripts/01_STEP6_merge_variants2.sh ${MASTER_DIR} ${MAPQUAL} ${CPUs}
        echo "log files saved in ${MASTER_DIR}/logs/bcftools_merge_logsIII"
fi

#STEP 7 - filter variants
if [[ "${step}" == 7 ]]; then
        echo "submitting the following job:"
        echo "sbatch ./01_subscripts/01_STEP7_filter_merged_variants2.sh ${MASTER_DIR} ${CPUs} ${CLEANUP} ${ID} ${MIND} ${MAF} ${QUAL} ${MIN_DEPTH} ${MAX_DEPTH} ${MAPQUAL}"
        mkdir -p ${MASTER_DIR}/logs/vcftools_filter_logsII
        mkdir -p ${MASTER_DIR}/filtered_biallelic_snpsII
        sbatch ./01_subscripts/01_STEP7_filter_merged_variants2.sh ${MASTER_DIR} ${CPUs} ${CLEANUP} ${ID} ${MIND} ${MAF} ${QUAL} ${MIN_DEPTH} ${MAX_DEPTH} ${MAPQUAL}
        echo "log files saved in ${MASTER_DIR}/logs/vcftools_filter_logsII"
fi




############################################################
# OUTPUTS                                                  #
############################################################
#TBD
