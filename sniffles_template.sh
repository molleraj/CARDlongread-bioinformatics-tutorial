#!/bin/bash

#SBATCH --cpus-per-task=16
#SBATCH --mem=16g
#SBATCH --mail-type=BEGIN,TIME_LIMIT_90,END
#SBATCH --time=96:00:00
#SBATCH --gres=lscratch:20

N=${SLURM_ARRAY_TASK_ID}
SAMPLE_SHEET='/data/CARD_AUX/LRS_temp/CHILE/SCRIPTS/CHILE_SAMPLE_FLOWCELL.tsv'
SAMPLE_ID=$(cut -f1 $SAMPLE_SHEET | sort | uniq | sed -n ${N}p)
# FLOWCELL=$(sed -n ${N}p $SAMPLE_SHEET | cut -f 2)

echo "Job Array #${N}"
echo "SAMPLE_ID ${SAMPLE_ID}"
# echo "FLOWCELL ${FLOWCELL}"

ml sniffles/2.3.3
BASE_DIR=/data/CARD_AUX/LRS_temp/CHILE

mkdir -p ${BASE_DIR}/Sniffles/snf

BAM=${BASE_DIR}/MAPPED_BAM/${SAMPLE_ID}.sorted_meth.bam
echo $BAM
sniffles --allow-overwrite -t 16 --input $BAM --snf ${BASE_DIR}/Sniffles/snf/${SAMPLE_ID}.snf

# combine VCFs
# ls /data/CARDPB/data/JAPAN/Sniffles/snf/* > SNIFFLES_samples.tsv
# mkdir -p /data/CARDPB/data/JAPAN/Sniffles/vcf
# sniffles --input SNIFFLES_samples.tsv --vcf JAPAN_PRKN_het.vcf

# split merged VCF by sample
# bcftools +split CHILE_sniffles_het.vcf -o Chile_sniffles_het_sample_split
