#!/usr/bin/bash
#SBATCH --cpus-per-task=40
#SBATCH --mem=120g
#SBATCH --mail-type=BEGIN,TIME_LIMIT_90,END
#SBATCH --time=24:00:00

#NEED TO CHANGE HERE
BASE_DIR=/data/CARDPB/data/RUSH
DEST_DIR=/data/CARD_AUX/LRS_temp/RUSH
SAMPLE_SHEET='/data/CARDPB/data/RUSH/SCRIPTS/RUSH_SAMPLE_FLOWCELL_UNIQUE_PAIRS_112724.tsv'

# add merged BAM directory
# mkdir -p ${BASE_DIR}/ONT_UBAM/

N=${SLURM_ARRAY_TASK_ID}

SAMPLE_ID=$(sed -n ${N}p $SAMPLE_SHEET | cut -f 1)
FLOWCELL=$(sed -n ${N}p $SAMPLE_SHEET | cut -f 2)


echo "Job Array #${N}"
echo "SAMPLE_ID ${SAMPLE_ID}"
echo "FLOWCELL ${FLOWCELL}"


echo "map with minimap2"
module load minimap2
module load samtools

mkdir -p ${DEST_DIR}/MAPPED_BAM/${SAMPLE_ID}

# index input unmapped BAM
samtools index ${BASE_DIR}/ONT_UBAM/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}_pod5.5mC_5hmC_v5.0.0.bam

# do mapping
samtools fastq -TMM,ML ${BASE_DIR}/ONT_UBAM/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}_pod5.5mC_5hmC_v5.0.0.bam | minimap2 -y -x map-ont --MD -t 20 -a --eqx -k 17 -K 10g /data/CARDPB/resources/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa - | samtools view -@ 10 -bh - | samtools sort -@ 10 - > ${DEST_DIR}/MAPPED_BAM/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}.sorted_meth.bam

# index mapped bam
samtools index ${DEST_DIR}/MAPPED_BAM/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}.sorted_meth.bam
