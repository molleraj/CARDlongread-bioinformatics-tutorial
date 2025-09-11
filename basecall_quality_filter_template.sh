#!/usr/bin/bash
#SBATCH --cpus-per-task=40
#SBATCH --mem=120g
#SBATCH --mail-type=BEGIN,TIME_LIMIT_90,END
#SBATCH --time=4:00:00

#NEED TO CHANGE HERE
BASE_DIR=/data/CARDPB/data/cohort
#DEST_DIR=/data/CARD_AUX/LRS_temp/UAB
SAMPLE_SHEET='/path/to/samplesheet'

N=${SLURM_ARRAY_TASK_ID}

SAMPLE_ID=$(sed -n ${N}p $SAMPLE_SHEET | cut -f 1)

# debugging output for logs
echo "Job Array #${N}"
echo "SAMPLE_ID ${SAMPLE_ID}"

# load modules
module load samtools nanopack minimap2

# make output directory if necessary
mkdir -p ${BASE_DIR}/MAPPED_BAM_BQ_FILTERED

# first convert UBAM to fastq, then filter by Q score, then map to hg38, then convert to bam, then sort
samtools fastq -TMm,Ml,MM,ML -@ 40 ${BASE_DIR}/ONT_UBAM/${SAMPLE_ID}_pod5.5mC_5hmC.bam | \
chopper -t 40 -q 10 | \
samtools import - -O bam -@ 10 | \
samtools sort -@ 10 > ${BASE_DIR}/MAPPED_BAM_BQ_FILTERED/${SAMPLE_ID}.sorted_bq10_filtered.bam

# make output mapped BAM index
samtools index ${BASE_DIR}/MAPPED_BAM_BQ_FILTERED/${SAMPLE_ID}.sorted_bq10_filtered.bam

# run cramino
BAM=${SAMPLE_ID}.sorted_bq10_filtered.bam

cramino -t 32 --reference /data/CARDPB/resources/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa ${BAM} > ${BAM}_cramino_QC.txt
