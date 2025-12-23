#!/usr/bin/env bash

#SBATCH --partition=gpu
#SBATCH --cpus-per-task=30
#SBATCH --mem=120g
#SBATCH --mail-type=BEGIN,TIME_LIMIT_90,END
#SBATCH --time=1:00:00
#SBATCH --gres=lscratch:20,gpu:a100:4

# get array job number for spooling subjobs
N=${SLURM_ARRAY_TASK_ID}
# specify base dir
BASE_DIR=/data/CARDPB/data/cohort
# specify path to sample sheet
# text list of pod5 group folders per sample/flowcell
# ls -d /data/CARDPB/data/PPMI/MethSmoothEval/POD5/*/*/channel_subset/*_dir
SAMPLE_SHEET='/path/to/samplesheet'
# first column is sample id (e.g., RUSH_001_FTX)
SAMPLE_ID=$(sed -n ${N}p $SAMPLE_SHEET | rev | cut -f4 -d'/' | rev)
# second column is flow cell (e.g., PAY78456)
FLOWCELL=$(sed -n ${N}p $SAMPLE_SHEET | rev | cut -f3 -d'/' | rev)
# third column is chunk (e.g., channel_subset_pod5_group_00)
CHUNK=$(sed -n ${N}p $SAMPLE_SHEET | rev | cut -f1 -d'/' | rev | sed 's/_dir//g')

# debugging output to slurm script
echo "Job Array #${N}"
echo "SAMPLE_ID ${SAMPLE_ID}"
echo "FLOWCELL ${FLOWCELL}"

# make output directory and parent if necessary
mkdir -p ${BASE_DIR}/ONT_UBAM/${SAMPLE_ID}

# load modules
# change 0.8.1 to 1.1.1 module (not yet default)
# change 1.1.1 to 0.9.6
# this was an example for 4 kHz basecalling that included a 4 kHz compatible model - change model to one desired for your use
module load dorado/0.9.6
module load pod5

# Basecall
# debugging output with output path for unmapped BAM
echo "${BASE_DIR}/ONT_UBAM/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}_${CHUNK}_pod5.5mC_5hmC_v4.1.0.bam"
dorado basecaller -x cuda:all ${DORADO_MODELS}/dna_r10.4.1_e8.2_400bps_sup@v4.1.0 ${BASE_DIR}/POD5/${SAMPLE_ID}/${FLOWCELL}/channel_subset/${CHUNK}_dir --skip-model-compatibility-check --modified-bases 5mCG_5hmCG > \
${BASE_DIR}/ONT_UBAM/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}_${CHUNK}_pod5.5mC_5hmC_v4.1.0.bam

# Summary
# debugging output with output path for sequencing summary QC text file
echo "${BASE_DIR}/ONT_UBAM/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}_${CHUNK}_sequencing_summary_v4.1.0.txt"
dorado summary ${BASE_DIR}/ONT_UBAM/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}_${CHUNK}_pod5.5mC_5hmC_v4.1.0.bam > ${BASE_DIR}/ONT_UBAM/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}_${CHUNK}_sequencing_summary_v4.1.0.txt
