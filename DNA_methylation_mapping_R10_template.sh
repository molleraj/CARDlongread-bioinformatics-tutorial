#!/usr/bin/bash
#SBATCH --cpus-per-task=40
#SBATCH --mem=120g
#SBATCH --mail-type=BEGIN,TIME_LIMIT_90,END
#SBATCH --time=24:00:00

#NEED TO CHANGE HERE
# base directory
BASE_DIR=/data/CARDPB/data/cohort
# in this case, also provide destination directory
DEST_DIR=/data/CARD_AUX/LRS_temp/cohort
# tab delimited file with list of samples and flow cells in columns 1 and 2
SAMPLE_SHEET='/path/to/samplesheet'

# add merged BAM directory
# mkdir -p ${BASE_DIR}/ONT_UBAM/

# get array job number for spooling subjobs
N=${SLURM_ARRAY_TASK_ID}

# first column is sample id (e.g., RUSH_001_FTX)
SAMPLE_ID=$(sed -n ${N}p $SAMPLE_SHEET | cut -f 1)
# second column is flow cell (e.g., PAY78456)
FLOWCELL=$(sed -n ${N}p $SAMPLE_SHEET | cut -f 2)

# debugging output to slurm script
echo "Job Array #${N}"
echo "SAMPLE_ID ${SAMPLE_ID}"
echo "FLOWCELL ${FLOWCELL}"

# debugging output
echo "map with minimap2"
# load modules (newest versions)
module load minimap2
module load samtools

# make output directory and parent if necessary
mkdir -p ${DEST_DIR}/MAPPED_BAM/${SAMPLE_ID}

# index input unmapped BAM
samtools index ${BASE_DIR}/ONT_UBAM/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}_pod5.5mC_5hmC_v5.0.0.bam

# do mapping
samtools fastq -TMM,ML ${BASE_DIR}/ONT_UBAM/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}_pod5.5mC_5hmC_v5.0.0.bam | minimap2 -y -x map-ont --MD -t 20 -a --eqx -k 17 -K 10g /data/CARDPB/resources/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa - | samtools view -@ 10 -bh - | samtools sort -@ 10 - > ${DEST_DIR}/MAPPED_BAM/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}.sorted_meth.bam

# index mapped bam
samtools index ${DEST_DIR}/MAPPED_BAM/${SAMPLE_ID}/${SAMPLE_ID}_${FLOWCELL}.sorted_meth.bam
