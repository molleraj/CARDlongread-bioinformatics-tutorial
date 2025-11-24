#!/bin/bash
# script to liftover microarray calls to reference of interest using dbSNP
# script depends on bcftools and plink; load as Biowulf modules
module load bcftools plink

# prompt for options
helpFunction()
{
   echo ""
   echo "Usage: $0 -d [dbSNP VCF] -r [reference name, e.g., hg38] -f [reference FASTA for reference allele correction] -i [input microarray Plink bed file (.bed/.bim/.fam)] -o [output lifted over microarray variants (.vcf.gz)]"
   echo -e "\t-d Path to dbSNP VCF for reference of interest"
   echo -e "\t-f Path to target reference for reference allele correction"
   echo -e "\t-r Target reference name (optional)"
   echo -e "\t-i Input microarray variant Plink bed file"
   echo -e "\t-o Output rsID-based liftover microarray variant file (.vcf.gz)"
   echo -e "\t-t Threads assigned to Plink and bcftools (optional)"
   exit 1 # Exit script after printing help
}

# get options listed above
while getopts "d:r:f:i:o:t:" opt
do
   case "$opt" in
      d ) dbsnp_variants="$OPTARG" ;;
      f ) reference_fasta="$OPTARG" ;;
      r ) reference_name="$OPTARG" ;;
      i ) input="$OPTARG" ;;
      o ) output="$OPTARG" ;;
      t ) threads="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case required parameters are empty
if [ -z "$dbsnp_variants" ] || [ -z "$reference_fasta"] || [ -z "$input" ] || [ -z "$output" ]
then
   echo "Some or all of the required parameters are empty.";
   helpFunction
fi

# Begin script in case all parameters are correct
# print out input variables for debugging
echo "dbSNP variants path: ${dbsnp_variants}"
echo "Reference fasta: ${reference_fasta}"
if [ ${reference_name} ]
then
	echo "Target reference: ${reference_name}"
fi
echo "Input microarray variants: ${input}"
echo "Output liftover: ${output}"
if [ ${threads} ]
then
	echo "Threads: ${threads}"
fi 

# get list of microarray rsIDs - always second column of Plink BED bim file
cut -f2 ${input}.bim > ${input}_microarray_variant_IDs_for_query.txt
# make copy and remove GSA- prefix in case array is GSA/GSAv3
cp ${input}_microarray_variant_IDs_for_query.txt ${input}_microarray_variant_IDs_corrected_names_for_query.txt
sed -i 's/GSA-//g' ${input}_microarray_variant_IDs_corrected_names_for_query.txt 
# save an ID correction file for Plink for later...
paste ${input}_microarray_variant_IDs_for_query.txt ${input}_microarray_variant_IDs_corrected_names_for_query.txt > ${input}_prefix_correction_table.tsv

# query hg38 dbSNP to get table with new coordinates for liftover
# probably rate limiting step in liftover - pretty slow to query all dbsnp
bcftools query -i ID==@${input}_microarray_variant_IDs_corrected_names_for_query.txt -f "%CHROM\t%POS\t%ID\t%REF\t%ALT" ${dbsnp_variants} > \
${input}_matches_corrected_names_to_dbSNP.tsv

# get rsID matches for later filtering in final liftover (only include rsIDs common between dbSNP and microarray batch1 for ACACB in final liftover)
cut -f3 ${input}_matches_corrected_names_to_dbSNP.tsv > ${input}_matches_corrected_names_to_dbSNP_rsIDs.tsv

# remember to convert GSA names back if necessary
# make sure to update id names, remove duplicates, and sort variants
# need to output in pfile (pgen + pvar) format to sort variants
plink2 --threads $--bfile batch1 --update-name ${input}_prefix_correction_table.tsv --rm-dup force-first --sort-vars --make-pgen --out ${input}_updated_names_unique_sorted_vars

# above two methods failed so try relabeling with plink - update chromosome and position first
plink2 --pfile ${input}_updated_names_unique_sorted_vars --extract ${input}_matches_corrected_names_to_dbSNP_rsIDs.tsv \
--update-chr ${input}_matches_corrected_names_to_dbSNP.tsv 1 3 \
--update-map ${input}_matches_corrected_names_to_dbSNP.tsv 2 3 --sort-vars --make-pgen \
--out ${input}_updated_names_matching_ids_dbSNP_liftover_coords_corrected

# second step with fa correction and removing bad snp bases
# make sure chr1...chrM naming scheme to match hg38 DV/clair3 calls (not 1..M)
plink2 --pfile ${input}_updated_names_matching_ids_dbSNP_liftover_coords_corrected \
--ref-from-fa ${reference_fasta} \
--snps-only just-acgt --export vcf bgz --output-chr chrM \
--out ${input}_updated_names_matching_ids_dbSNP_liftover

# index output
bcftools index ${input}_updated_names_matching_ids_dbSNP_liftover.vcf.gz

# fix alleles relative to reference
bcftools +fixref ${input}_updated_names_matching_ids_dbSNP_liftover.vcf.gz -O z -o ${input}_updated_names_matching_ids_dbSNP_liftover_reffixed.vcf.gz \
-- -f ${reference_fasta} -m flip -d
# index output
bcftools index ${input}_updated_names_matching_ids_dbSNP_liftover_reffixed.vcf.gz
