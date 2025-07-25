#!/bin/bash

# Load parameters
DIR="$1"
ANNO="$2"
TXFILE="$3"
ANALYSIS_DIR="$4"
CLAIR_MODEL="$5"
MIN_Q="$6"
MAX_U="$7"
MAPQ="$8"
EXT="$9"
BED="${10}"

#pre-processing
mapfile -t SAMPLES < <(cut -f1 $ANNO)

QUANT=$((100 - MAX_U))
MOD="_q${QUANT}_Q${MIN_Q}"
MOD_STRP="q${QUANT}_Q${MIN_Q}"
MAPQ_MOD="_MAPQ${MAPQ}"
# Activate the required conda environment


OUT_DIR="$ANALYSIS_DIR/output/"
mkdir -p "$OUT_DIR"
mkdir -p "$OUT_DIRdepth/"

#######################


echo "################# BAM DEPTH PLOTTING ##################"
echo "This might take a while."


# plot depth
ARGS=("$ANALYSIS_DIR/filtered_bam_sr/depth/" "$EXT" "$MOD$MAPQ_MOD" "$ANNO" "$BED" "${OUT_DIR}depth/")
Rscript "./scripts/depth_analysis.R" "${ARGS[@]}"

echo "Per sample read depth plots and a summary of the depths statistics saved in: ${OUT_DIR}depth/. More detailed data lies in $ANALYSIS_DIR/filtered_bam_sr/depth/."

###########################


# run for snp and indel individually
Rscript ./scripts/vcf_amplicon_geno_clairs-to.R "$ANALYSIS_DIR/ClairS-TO/" "$MOD_STRP${MAPQ_MOD}_$CLAIR_MODEL" "$TXFILE" "$OUT_DIR" "snv" "${SAMPLES[@]}"
Rscript ./scripts/vcf_amplicon_geno_clairs-to.R "$ANALYSIS_DIR/ClairS-TO/" "$MOD_STRP${MAPQ_MOD}_$CLAIR_MODEL" "$TXFILE" "$OUT_DIR" "indel" "${SAMPLES[@]}"


echo "SNV genotyping info saved in: $ANALYSIS_DIR/ClairS-TO/$MOD_STRP${MAPQ_MOD}_${CLAIR_MODEL}_snv_vcf_collection.tsv."
echo "Indel genotyping info saved in: $ANALYSIS_DIR/ClairS-TO/$MOD_STRP${MAPQ_MOD}_${CLAIR_MODEL}_indel_vcf_collection.tsv."

Rscript ./scripts/genotyping.R "$ANALYSIS_DIR/ClairS-TO/" "${MOD_STRP}${MAPQ_MOD}_${CLAIR_MODEL}_snv_vcf_collection.tsv" "$OUT_DIR"

echo "Genotyping results saved in: ${OUT_DIR}SNV_genotyping_results.tsv, ${OUT_DIR}Prot_coding_SNV_genotyping_results.tsv and ${OUT_DIR}Indel_genotyping_results.tsv"

exit 0
