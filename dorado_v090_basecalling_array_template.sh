#!/usr/bin/env bash

#SBATCH --partition=gpu
#SBATCH --cpus-per-task=30
#SBATCH --mem=120g
#SBATCH --mail-type=BEGIN,TIME_LIMIT_90,END
#SBATCH --time=48:00:00
#SBATCH --gres=lscratch:50,gpu:a100:2

# get array job number for spooling subjobs
N=${SLURM_ARRAY_TASK_ID}
# specify path to sample sheet
# tab delimited file with list of samples and flow cells in columns 1 and 2
SAMPLE_SHEET='/path/to/samplesheet'
# first column is sample id (e.g., RUSH_001_FTX)
SAMPLE_ID=$(sed -n ${N}p $SAMPLE_SHEET | cut -f 1)
# second column is flow cell (e.g., PAY78456)
FLOWCELL=$(sed -n ${N}p $SAMPLE_SHEET | cut -f 2)

# debugging output to slurm script
echo "Job Array #${N}"
echo "SAMPLE_ID ${SAMPLE_ID}"
echo "FLOWCELL ${FLOWCELL}"

# base directory path
BASE_DIR=/data/CARDPB/data/cohort # for example

# make output directory and parent if necessary
mkdir -p ${BASE_DIR}/ONT_UBAM/${SAMPLE_ID}

# load modules
# change 0.8.1 to 0.9.0 module (not yet default)
module load dorado/0.9.0
module load pod5

# Basecall
# debugging output with output path for unmapped BAM
echo "${BASE_DIR}/ONT_UBAM/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}_pod5.5mC_5hmC_v5.0.0.bam"
dorado basecaller -x cuda:all ${DORADO_MODELS}/dna_r10.4.1_e8.2_400bps_sup@v5.0.0 ${BASE_DIR}/POD5/${SAMPLE_ID}/${FLOWCELL} --skip-model-compatibility-check --modified-bases 5mC_5hmC > ${BASE_DIR}/ONT_UBAM/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}_pod5.5mC_5hmC_v5.0.0.bam

# Summary
# debugging output with output path for sequencing summary QC text file
echo "${BASE_DIR}/ONT_UBAM/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}_sequencing_summary_v5.0.0.txt"
dorado summary ${BASE_DIR}/ONT_UBAM/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}_pod5.5mC_5hmC_v5.0.0.bam > ${BASE_DIR}/ONT_UBAM/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}_sequencing_summary_v5.0.0.txt
