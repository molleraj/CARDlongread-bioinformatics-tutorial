#!/bin/bash

#SBATCH --cpus-per-task=16
#SBATCH --mem=32g
#SBATCH --mail-type=BEGIN,TIME_LIMIT_90,END
#SBATCH --time=1:00:00
#SBATCH --gres=lscratch:20

# get array job number for spooling subjobs
N=${SLURM_ARRAY_TASK_ID}
# set base directory
BASE_DIR=/data/CARDPB/data/cohort
# tab delimited file with list of samples and flow cells in columns 1 and 2
SAMPLE_SHEET='/path/to/samplesheet'
SAMPLE_ID=$(sed -n ${N}p $SAMPLE_SHEET| cut -f 1)
FLOWCELL=$(sed -n ${N}p $SAMPLE_SHEET | cut -f 2)

# debugging output to slurm script
echo "Job Array #${N}"
echo "SAMPLE_ID ${SAMPLE_ID}"
echo "FLOWCELL ${FLOWCELL}"

# load annovar for SNV annotation
ml annovar

# set VCF variable for running AnnoVar based on per sample concatenated DV VCF

VCF=${BASE_DIR}/DEEPVARIANT/${SAMPLE_ID}/${SAMPLE_ID}_deepvariant_concat.vcf.gz

# annotation based on Pilar's previously documented script
# part 2 of https://github.com/neurogenetics/NBIA_PD

# Annotating variants
table_annovar.pl ${VCF} $ANNOVAR_DATA/hg38 -buildver hg38 \
--thread 16 \
-out ${BASE_DIR}/ANNOVAR/${SAMPLE_ID}.annovar \
-remove -protocol ensGene,ljb26_all,gnomad211_genome,clinvar_20190305 \
-operation g,f,f,f \
-nastring . \
-vcfinput
