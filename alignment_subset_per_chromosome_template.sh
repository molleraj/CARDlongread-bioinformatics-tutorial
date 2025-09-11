#!/usr/bin/env bash

#SBATCH --partition=norm
#SBATCH --cpus-per-task=64
#SBATCH --mem=64g
#SBATCH --time=6:00:00
#SBATCH --gres=lscratch:10
#SBATCH --mail-type=BEGIN,TIME_LIMIT_90,END

#NEED TO CHANGE HERE
BASE_DIR=/data/CARDPB/data/E46K
#DEST_DIR=/data/CARD_AUX/LRS_temp/RUSH
SAMPLE_SHEET='/data/CARDPB/data/E46K/SCRIPTS/E46K_SAMPLE_FLOWCELL.tsv'

# add merged BAM directory
# mkdir -p ${BASE_DIR}/ONT_UBAM_CARDPB/

N=${SLURM_ARRAY_TASK_ID}

SAMPLE_ID=$(sed -n ${N}p $SAMPLE_SHEET | cut -f 1)
FLOWCELL=$(sed -n ${N}p $SAMPLE_SHEET | cut -f 2)

echo "Job Array #${N}"
echo "SAMPLE_ID ${SAMPLE_ID}"
# echo "FLOWCELL ${FLOWCELL}"

# source myconda
# conda activate base

module load samtools

# chromosome list
CHROMOSOMES=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY" "chrM")

# specify mapped BAM input
BAM_IN=${BASE_DIR}/MAPPED_BAM/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}.second_attempt.sorted_meth.bam
for i in "${CHROMOSOMES[@]}"
do
BAM_SUBSET=${SAMPLE_ID}_${FLOWCELL}.${i}.second_attempt.sorted_meth.bam
# specify output directory as SAMPLE ID parent and FLOWCELL descendant
OUTPUT=${BASE_DIR}/MAPPED_BAM/${SAMPLE_ID}_per_chromosome
# DIR_NAME=$(dirname ${BAM_IN})
# specify sample name for VCF file
# SAMPLENAME=${SAMPLE_ID}_${FLOWCELL}

mkdir -p ${OUTPUT}
# make directory and parent accessible to everyone on CARDPB
# chmod -R 775 ${BASE_DIR}/CLAIR3_CHR20/${SAMPLE_ID}
chmod -R 775 ${OUTPUT}

# make subset BAM
samtools view -b --threads 60 ${BAM_IN} ${i} > ${OUTPUT}/${BAM_SUBSET}
# index BAM subset
samtools index ${OUTPUT}/${BAM_SUBSET}

# clair3 --bam_fn=${BAM_SUBSET} --ref_fn=/data/CARDPB/resources/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa --threads=24 --platform=ont --model_path=${CLAIR3_MODELS}/ont/ --output=${OUTPUT} --sample_name=${SAMPLENAME}

# make directory AND output files accessible to everyone on CARDPB
chmod -R 775 ${OUTPUT}
done

# # how to run
# # sbatch --mem=10g --cpus-per-task=2 --time=4-0 clair3.sh /data/CARDPB/data/HBCC/MAPPED_HG38/"$line"_HG38.sorted.bam /data/CARDPB/data/HBCC/VARIANTS/CLAIR3/PER_FC/"$line
