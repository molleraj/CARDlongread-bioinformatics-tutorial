#!/usr/bin/env bash

#SBATCH --partition=norm
#SBATCH --cpus-per-task=30
#SBATCH --mem=30g
#SBATCH --mail-type=BEGIN,TIME_LIMIT_90,END
#SBATCH --time=2:00:00
#SBATCH --gres=lscratch:20

# get array job number for spooling subjobs
N=${SLURM_ARRAY_TASK_ID}
# specify base dir
BASE_DIR=/data/CARDPB/data/cohort
# specify path to sample sheet
# text list of pod5 group folders per sample/flowcell
# ls -d /data/CARDPB/data/MOSR/POD5/*/*/channel_subset/*_dir
SAMPLE_SHEET='/data/CARDPB/data/cohort/SCRIPTS/cohort_sample_flowcell.tsv'
# first column is sample id (e.g., RUSH_001_FTX)
SAMPLE_ID=$(sed -n ${N}p $SAMPLE_SHEET | cut -f1)
# second column is flow cell (e.g., PAY78456)
FLOWCELL=$(sed -n ${N}p $SAMPLE_SHEET | cut -f2)

# debugging output to slurm script
echo "Job Array #${N}"
echo "SAMPLE_ID ${SAMPLE_ID}"
echo "FLOWCELL ${FLOWCELL}"

# make output directory and parent if necessary
mkdir -p ${BASE_DIR}/ONT_UBAM/${SAMPLE_ID}

# load samtools module
module load samtools

# this time use sample integrity check with samtools quickcheck
ls ${BASE_DIR}/ONT_UBAM/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}_channel_subset_pod5_group*_pod5.5mCG_5hmCG_v5.2.0.bam | \
	grep -v -f <(samtools quickcheck -u -v ${BASE_DIR}/ONT_UBAM/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}_channel_subset_pod5_group*_pod5.5mCG_5hmCG_v5.2.0.bam) > \
	${BASE_DIR}/ONT_UBAM/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}_successful_channel_subset_pod5_ubams.txt
# merge channel groups into single BAM
samtools merge -f -@ 60 -o ${BASE_DIR}/ONT_UBAM/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}_merged_pod5.5mCG_5hmCG_v5.2.0.bam \
	-b ${BASE_DIR}/ONT_UBAM/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}_successful_channel_subset_pod5_ubams.txt
# index BAM for downstream steps
samtools index -@ 60 ${BASE_DIR}/ONT_UBAM/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}_merged_pod5.5mCG_5hmCG_v5.2.0.bam
