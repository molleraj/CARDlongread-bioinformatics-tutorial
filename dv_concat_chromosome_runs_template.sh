#!/usr/bin/env bash
# script for preparing SRS VCF subsets from AMP-PD VCFs Ken provided
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=30g
#SBATCH --mail-type=BEGIN,TIME_LIMIT_90,END
#SBATCH --time=1:00:00
# bcftools view --threads 30 -S /data/CARDPB/data/PPMI/scripts/PPMI_chr20_amp_pd_VCF_sample_list.txt chr20.vcf.gz -O z -o /data/CARDPB/data/PPMI/SRS_PPMI_VCF_SUBSETS/PPMI_chr20.vcf.gz
# N=${SLURM_ARRAY_TASK_ID}
# VCF_LIST=/data/CARDPB/data/PPMI/scripts/AMP-PD_VCF_list_102124.txt
# AMP_PD_VCF_PATH=/data/CARD/PD/AMP-PD/VCFs
# CHR_VCF=$(sed -n ${N}p $VCF_LIST)

module load samtools bcftools

# get array job ID
N=${SLURM_ARRAY_TASK_ID}
# echo "Job Array #${N}"
# echo "CHR_VCF ${CHR_VCF}"

# bcftools view --threads 60 -S /data/CARDPB/data/PPMI/scripts/PPMI_chr20_amp_pd_VCF_sample_list.txt ${AMP_PD_VCF_PATH}/${CHR_VCF} -O z -o /data/CARDPB/data/PPMI/SRS_PPMI_VCF_SUBSETS/PPMI_${CHR_VCF}
# concatenate RUSH vcfs in /data/CARD_AUX
BASE_DIR=/data/CARDPB/data/cohort
# SRS_VCFS_PATH=/data/CARDPB/data/PPMI/SRS_PPMI_VCF_SUBSETS
SAMPLE_SHEET='/path/to/samplesheet'
# get sample ID from sample sheet
SAMPLE_ID=$(sed -n ${N}p $SAMPLE_SHEET | cut -f1)
FLOWCELL=$(sed -n ${N}p $SAMPLE_SHEET | cut -f2)

# debugging output
echo "Job Array #${N}"
echo "Sample ${SAMPLE_ID}"
echo "Flow cell ${FLOWCELL}"

# concat chromosome 20 vcf with clair3 chromosome 20
# bcftools index --threads 60 PPMI_chr20_clair3_vcf_merged_101724.vcf.gz
# bcftools index --threads 60 ${SRS_VCFS_PATH}/PPMI_chr20.vcf.gz
# bcftools merge --threads 60 ${SRS_VCFS_PATH}/PPMI_chr20.vcf.gz PPMI_chr20_clair3_vcf_merged_101724.vcf.gz -O z -o ${SRS_VCFS_PATH}/PPMI_chr20_SRS_LRS_merge_102124.vcf.gz
# concat WGS vcf with deepvariant WGS
bcftools concat --threads 30 ${BASE_DIR}/DEEPVARIANT/${SAMPLE_ID}/*/*FB.vcf.gz -O z -o ${BASE_DIR}/DEEPVARIANT/${SAMPLE_ID}/${SAMPLE_ID}_deepvariant_concat.vcf.gz
bcftools index ${BASE_DIR}/DEEPVARIANT/${SAMPLE_ID}/${SAMPLE_ID}_deepvariant_concat.vcf.gz
# correctly rename samples
# bcftools reheader --threads 60 -s PPMI_SRS_LRS_corrected_sample_names_for_VCF_102524.txt ${SRS_VCFS_PATH}/PPMI_WGS_SRS_LRS_merge_102124.vcf.gz -o ${SRS_VCFS_PATH}/PPMI_WGS_SRS_LRS_merge_corrected_names_102524.vcf.gz
