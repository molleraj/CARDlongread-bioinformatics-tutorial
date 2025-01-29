#!/bin/bash

#SBATCH --cpus-per-task=16
#SBATCH --mem=16g
#SBATCH --mail-type=BEGIN,TIME_LIMIT_90,END
#SBATCH --time=1:00:00
#SBATCH --gres=lscratch:20

# get array job number for spooling subjobs
N=${SLURM_ARRAY_TASK_ID}
# tab delimited file with list of samples and flow cells in columns 1 and 2
SAMPLE_SHEET='/path/to/samplesheet'
SAMPLE_ID=$(sed -n ${N}p $SAMPLE_SHEET| cut -f 1)
# FLOWCELL=$(sed -n ${N}p $SAMPLE_SHEET | cut -f 2)

# debugging output to slurm script
echo "Job Array #${N}"
echo "SAMPLE_ID ${SAMPLE_ID}"
# echo "FLOWCELL ${FLOWCELL}"

# load sniffles version 2.5.3
ml sniffles/2.5.3
BASE_DIR=/data/CARDPB/data/cohort

# make snf output directory and parents if necessary
mkdir -p ${BASE_DIR}/Sniffles/snf

BAM=${BASE_DIR}/MAPPED_BAM/${SAMPLE_ID}.sorted_meth.bam
echo $BAM
# run sniffles with 16 threads 
sniffles --allow-overwrite -t 16 --input $BAM --snf ${BASE_DIR}/Sniffles/snf/${SAMPLE_ID}.snf

# combine VCFs
ls ${BASE_DIR}/Sniffles/snf/* > SNIFFLES_samples.tsv
mkdir -p ${BASE_DIR}/Sniffles/vcf
OUTPUT_VCF=cohort_sniffles_het.vcf
sniffles --input SNIFFLES_samples.tsv --vcf ${OUTPUT_VCF}

# split merged VCF by sample
ml samtools
bcftools +split cohort_sniffles_het.vcf -o cohort_het_sample_split
