#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=32g
#SBATCH --time=6:00:00

N=${SLURM_ARRAY_TASK_ID}
SAMPLE_SHEET='/data/CARD_AUX/LRS_temp/CHILE/SCRIPTS/CHILE_SAMPLE_FLOWCELL.tsv'
SAMPLE_ID=$(cut -f1 $SAMPLE_SHEET | sort | uniq | sed -n ${N}p)
# FLOWCELL=$(sed -n ${N}p $SAMPLE_SHEET | cut -f 2)

echo "Job Array #${N}"
echo "SAMPLE_ID ${SAMPLE_ID}"
# echo "FLOWCELL ${FLOWCELL}"

BASE_DIR=/data/CARD_AUX/LRS_temp/CHILE

BAM=${BASE_DIR}/MAPPED_BAM/${SAMPLE_ID}.sorted_meth.bam
OUTPUT_VCF=${BASE_DIR}/VAMOS/${SAMPLE_ID}_vamos.vcf
SAMPLENAME=${SAMPLE_ID}

# Run with docker

VAMOS=/data/CARDPB/data/adaptive_sampling/Madison_repeat/VAMOS/vamos_2.1.3.sif

ml singularity

mkdir -p ${BASE_DIR}/VAMOS

#connect files (CHANGE HERE when change the datasets!)
export SINGULARITY_BINDPATH="/data/CARD_AUX/LRS_temp/CHILE/VAMOS,/data/CARD_AUX/LRS_temp/CHILE/MAPPED_BAM,/data/CARDPB/data/adaptive_sampling/Madison_repeat/VAMOS"

# /data/CARDPB/data/NABEC/repeat_expression/vntrs_motifs_delta_0.2.bed >> suggested to use for motif _analysis from Dr.Mark Chaisson.
# /data/CARD/FOUNDINtemp/LRS_QTLs/vamos_expression/references_vntr_annotations_q-0.1_motifs.tsv   >> in Terra NABEC samples/

singularity exec $VAMOS vamos --read -b ${BAM} \
-r /data/CARDPB/data/adaptive_sampling/Madison_repeat/VAMOS/vamos.motif.hg38.v2.1.e0.1.tsv \
-s ${SAMPLENAME} \
-o ${OUTPUT_VCF} \
-t 4

echo Done!!
