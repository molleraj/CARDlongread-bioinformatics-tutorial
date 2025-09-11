#!/bin/bash

#SBATCH --partition=norm
#SBATCH --cpus-per-task=16
#SBATCH --mem=64g
#SBATCH --mail-type=BEGIN,TIME_LIMIT_90,END
#SBATCH --time=4:00:00
#SBATCH --gres=lscratch:50

N=${SLURM_ARRAY_TASK_ID}
# note this has a different structure from a typical sample sheet
# list of subsetted BAM files
SAMPLE_SHEET='/path/to/samplesheet'
BAM_SUBSET=$(sed -n ${N}p $SAMPLE_SHEET)
SAMPLE_ID=$(sed -n ${N}p $SAMPLE_SHEET | cut -f1 -d'.')
CHROMOSOME=$(sed -n ${N}p $SAMPLE_SHEET | cut -f2 -d'.')

# debugging output
echo "Job Array #${N}"
echo "BAM_SUBSET ${BAM_SUBSET}"
echo "SAMPLE_ID ${SAMPLE_ID}"
echo "CHROMOSOME ${CHROMOSOME}"

# specify base directory
BASE_DIR=/data/CARDPB/data/cohort

# module load singularity
module load deepvariant/1.8.0

# make DEEPVARIANT output folder and parents if necessary
mkdir -p ${BASE_DIR}/DEEPVARIANT

# use per sample per chromosome subsetted BAMs
BAM_IN=${BASE_DIR}/MAPPED_BAM/${SAMPLE_ID}_per_chromosome/${BAM_SUBSET}
OUT_FOLDER=${BASE_DIR}/DEEPVARIANT/${SAMPLE_ID}/${CHROMOSOME}
SAMPLE_PREFIX=${SAMPLE_ID}
# deepvariant docker image version
# BIN_VERSION="1.8.0-gpu"

# make output directory and parent if necessary
mkdir -p ${OUT_FOLDER}
# cp -r $PEPPER_DATA/* .
# use local singularity/docker image
# docker://google/deepvariant:"${BIN_VERSION}-gpu" \
# change to local instance of docker image
# singularity run --nv -B /usr/lib/locale/:/usr/lib/locale/ \
# /data/mollerabg/deepvariant/deepvariant_1.8.0-gpu.sif \

# removed docker image described above. now use built-in biowulf module
run_deepvariant \
--model_type=ONT_R104 \
--ref=/data/CARDPB/resources/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
--reads=${BAM_IN} \
--sample_name=${SAMPLE_PREFIX} \
--output_vcf=${OUT_FOLDER}/${SAMPLE_PREFIX}.vcf.gz \
--output_gvcf=${OUT_FOLDER}/${SAMPLE_PREFIX}.g.vcf.gz \
--intermediate_results_dir=${OUT_FOLDER}/intermediate_results_dir \
--logging_dir=${OUT_FOLDER}/logs \
--call_variants_extra_args="writer_threads=1" \
--num_shards=16

# make output folder accessible to everybody
chmod -R 775 ${OUT_FOLDER}
# __________________________________________________________________________________________________
# INDEX PHASED BAM
# module load samtools
# samtools index ${OUT_FOLDER}/${SAMPLE_PREFIX}.haplotagged.bam
