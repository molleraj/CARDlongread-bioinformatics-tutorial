#!/usr/bin/env bash

#SBATCH --cpus-per-task=24
#SBATCH --mem=80g
#SBATCH --mail-type=BEGIN,TIME_LIMIT_90,END
#SBATCH --time=8:00:00

# load modkit module
ml modkit

#NEED TO CHANGE HERE
BASE_DIR=/data/CARDPB/data/cohort
SAMPLE_SHEET='/path/to/samplesheet'

# add methylation BED directory
# mkdir -p ${BASE_DIR}/METHYLATION/

N=${SLURM_ARRAY_TASK_ID}

# get sample name and flow cell ID variables from sample sheet

SAMPLE_ID=$(sed -n ${N}p $SAMPLE_SHEET | cut -f 1)
FLOWCELL=$(sed -n ${N}p $SAMPLE_SHEET | cut -f 2)

mkdir -p ${BASE_DIR}/METHYLATION/${SAMPLE_ID}_${FLOWCELL}

echo "Job Array #${N}"
echo "SAMPLE_ID ${SAMPLE_ID}"
echo "FLOWCELL ${FLOWCELL}"

BAM_FILE=${BASE_DIR}/PEPPER_phased_bam/${SAMPLE_ID}_${FLOWCELL}/${SAMPLE_ID}_${FLOWCELL}.haplotagged.bam
OUT=${BASE_DIR}/METHYLATION/${SAMPLE_ID}_${FLOWCELL}

# use modkit to get methylation BAM
modkit pileup --cpg --ref /data/CARDPB/resources/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa --only-tabs --threads 24 --ignore h --combine-strands --partition-tag HP --prefix ${SAMPLE_ID}_${FLOWCELL} ${BAM_FILE} ${OUT}

# --include-bed ${BED_file_w_regions}
