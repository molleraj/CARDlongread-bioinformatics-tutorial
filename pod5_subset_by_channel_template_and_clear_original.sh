#!/usr/bin/env bash

#SBATCH --partition=norm
#SBATCH --cpus-per-task=56
#SBATCH --mem=40g
#SBATCH --mail-type=BEGIN,TIME_LIMIT_90,END
#SBATCH --time=2:00:00

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
# mkdir -p ${BASE_DIR}/POD5/${SAMPLE_ID}/${FLOWCELL}/channel_subset/

# load conda environment
source /data/$USER/conda/etc/profile.d/conda.sh
conda activate pod5
# we use pod5 version 0.3.27 for MinKNOW 24+ generated pod5s
# set ulimits to enough processes can run simultaneously
ulimit -n 32768 -u 32768

# get summary file
pod5 view ${BASE_DIR}/POD5/${SAMPLE_ID}/${FLOWCELL}/ --include "read_id, channel" --output ${BASE_DIR}/POD5/${SAMPLE_ID}/${FLOWCELL}/read_id_channel_summary.tsv --force-overwrite

# subset using summary file
# remove original POD5s if and only if subset completes successfully
pod5 subset --threads 24 ${BASE_DIR}/POD5/${SAMPLE_ID}/${FLOWCELL}/ --summary ${BASE_DIR}/POD5/${SAMPLE_ID}/${FLOWCELL}/read_id_channel_summary.tsv \
--columns channel --output ${BASE_DIR}/POD5/${SAMPLE_ID}/${FLOWCELL}/channel_subset --force-overwrite && rm -rf ${BASE_DIR}/POD5/${SAMPLE_ID}/${FLOWCELL}/*.pod5

# move subsets into appropriate directories
# make file list
ls ${BASE_DIR}/POD5/${SAMPLE_ID}/${FLOWCELL}/channel_subset/*.pod5 > ${BASE_DIR}/POD5/${SAMPLE_ID}/${FLOWCELL}/channel_subset/channel_subset_pod5_file_list.txt

# split into x parts (try 32)
split -d -n l/32 ${BASE_DIR}/POD5/${SAMPLE_ID}/${FLOWCELL}/channel_subset/channel_subset_pod5_file_list.txt ${BASE_DIR}/POD5/${SAMPLE_ID}/${FLOWCELL}/channel_subset/channel_subset_pod5_group_

# now move into appropriate directories
for i in $(find ${BASE_DIR}/POD5/${SAMPLE_ID}/${FLOWCELL}/channel_subset/channel_subset_pod5_group_* -type f)
do
mkdir ${i}_dir
mv $(cat ${i}) ${i}_dir
done

# next move to highly parallel basecalling
