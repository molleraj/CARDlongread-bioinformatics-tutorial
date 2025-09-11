#!/bin/bash

#SBATCH --cpus-per-task=16
#SBATCH --mem=16g
#SBATCH --mail-type=BEGIN,TIME_LIMIT_90,END
#SBATCH --time=1:00:00
#SBATCH --gres=lscratch:20

# get array job number for spooling subjobs
N=${SLURM_ARRAY_TASK_ID}
# set base directory; cohort name example below
BASE_DIR=/data/CARDPB/data/cohort
# tab delimited file with list of samples and flow cells in columns 1 and 2
SAMPLE_SHEET='/path/to/samplesheet'
SAMPLE_ID=$(sed -n ${N}p $SAMPLE_SHEET| cut -f 1)
FLOWCELL=$(sed -n ${N}p $SAMPLE_SHEET | cut -f 2)

# debugging output to slurm script
echo "Job Array #${N}"
echo "SAMPLE_ID ${SAMPLE_ID}"
echo "FLOWCELL ${FLOWCELL}"

# load sniffles version 2.5.3
ml annotsv bedtools bcftools

# set VCF variable for running AnnotSV based on Sniffles split VCF
# BAM=${BASE_DIR}/MAPPED_BAM/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}.second_attempt.sorted_meth.bam
# echo $BAM

# split merged VCF by sample
# ml samtools
# bcftools +split ${BASE_DIR}/Sniffles/vcf/${OUTPUT_VCF} -o ${BASE_DIR}/Sniffles/vcf/cohort_het_sample_split
VCF=${BASE_DIR}/Sniffles/vcf/cohort_het_sample_split/${SAMPLE_ID}.vcf.gz
# run AnnotSV and output to E46K_het_sample_split directory
AnnotSV -SVinputFile ${VCF} -SVinputInfo 1 -outputFile ${BASE_DIR}/Sniffles/vcf/cohort_het_sample_split/${SAMPLE_ID}.annotsv.tsv
