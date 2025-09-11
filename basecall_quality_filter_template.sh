#!/usr/bin/bash
#SBATCH --cpus-per-task=40
#SBATCH --mem=120g
#SBATCH --mail-type=BEGIN,TIME_LIMIT_90,END
#SBATCH --time=24:00:00

#NEED TO CHANGE HERE
BASE_DIR=/data/CARDPB/data/COLOMBIA
#DEST_DIR=/data/CARD_AUX/LRS_temp/UAB
SAMPLE_SHEET='/data/CARDPB/data/COLOMBIA/SCRIPTS/Colombia_ubam_merged_list_072125.txt'

N=${SLURM_ARRAY_TASK_ID}

SAMPLE_ID=$(sed -n ${N}p $SAMPLE_SHEET | cut -f 1)

# debugging output for logs
echo "Job Array #${N}"
echo "SAMPLE_ID ${SAMPLE_ID}"

# load modules
module load samtools nanopack minimap2

# make output directory if necessary
mkdir -p ${BASE_DIR}/MAPPED_BAM_MERGED_BQ_FILTERED

# first convert UBAM to fastq, then filter by Q score, then map to hg38, then convert to bam, then sort
samtools fastq -TMm,Ml,MM,ML -@ 40 ${BASE_DIR}/UBAM_MERGED/${SAMPLE_ID}_pod5.5mCG_5hmCG.bam | \
chopper -t 40 -q 10 | \
minimap2 -ax map-ont /data/CARDPB/resources/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa - -k 17 -y -K 10g -t 20 --MD --eqx | \
samtools view -@ 10 -bh - | \
samtools sort -@ 10 > ${BASE_DIR}/MAPPED_BAM_MERGED_BQ_FILTERED/${SAMPLE_ID}.sorted_minimap2_bq10_filtered.bam

# make output mapped BAM index
samtools index ${BASE_DIR}/MAPPED_BAM_MERGED_BQ_FILTERED/${SAMPLE_ID}.sorted_minimap2_bq10_filtered.bam

# run cramino
BAM=${SAMPLE_ID}

cramino -t 32 --reference /data/CARDPB/resources/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa ${BAM} > ${BAM}_cramino_QC.txt
