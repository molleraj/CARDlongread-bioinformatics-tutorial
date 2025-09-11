#!/bin/bash
# request interactive node with 32 CPUs and 32GB RAM
# sinteractive --cpus-per-task=32 --mem=32g
# never mind - run as a batch job
#SBATCH --cpus-per-task=32
#SBATCH --mem=32g
#SBATCH --time=2:00:00
# load conda environment
source /data/$USER/conda/etc/profile.d/conda.sh
conda activate tryvamos
# set BASE_DIR
BASE_DIR=/data/CARDPB/data/cohort
# combine VCFs
python ~/vamos/tryvamos/tryvamos.py combineVCF vamos_sample_vcf_input_list.txt ${BASE_DIR}/VAMOS/vamos_merged_cohort_variants.vcf
