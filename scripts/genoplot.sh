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

#preprocessing
mapfile -t SAMPLES < <(cut -f1 $ANNO)

QUANT=$((100 - MAX_U))
MOD="_q${QUANT}_Q${MIN_Q}"
MOD_STRP="q${QUANT}_Q${MIN_Q}"
MAPQ_MOD="_MAPQ${MAPQ}"
# Activate the required Conda environment

#conda deactivate
#conda activate nagger_plotting

OUT_DIR="$ANALYSIS_DIR/output/"
mkdir -p "$OUT_DIR"

#######################


echo "################# BAM DEPTH PLOTTING ##################"
echo "This might take a while."


# run for post-filtering data
ARGS=("$ANALYSIS_DIR/filtered_bam_sr/depth/" "$EXT" "$MOD$MAPQ_MOD" "$ANNO" "$BED" "$OUT_DIR")
Rscript "./scripts/depth_analysis.R" "${ARGS[@]}"

echo "Post-filtering per sample read depth plots saved in: $DIR/filtered_bam_sr/depth/plots as post-filtering_depth.svg"

###########################


# run for snp and indel individually
Rscript ./scripts/vcf_amplicon_geno_clairs-to.R "$ANALYSIS_DIR/ClairS-TO/" "$MOD_STRP${MAPQ_MOD}_$CLAIR_MODEL" "$TXFILE" "$OUT_DIR" "snv" "${SAMPLES[@]}"
Rscript ./scripts/vcf_amplicon_geno_clairs-to.R "$ANALYSIS_DIR/ClairS-TO/" "$MOD_STRP${MAPQ_MOD}_$CLAIR_MODEL" "$TXFILE" "$OUT_DIR" "indel" "${SAMPLES[@]}"


echo "SNV genotyping info saved in: $OUT_DIR$MOD_STRP${MAPQ_MOD}_${CLAIR_MODEL}_snv_vcf_collection.tsv."
echo "Indel genotyping info saved in: $OUT_DIR$MOD_STRP${MAPQ_MOD}_${CLAIR_MODEL}_indel_vcf_collection.tsv."

Rscript ./scripts/genotyping.R "$ANALYSIS_DIR/ClairS-TO/" "${MOD_STRP}${MAPQ_MOD}_${CLAIR_MODEL}_snv_vcf_collection.tsv" "$MOD_STRP" "$OUT_DIR"

echo "Genotyping results saved in: $OUT_DIR$MOD_STRP${MAPQ_MOD}-genotyping_results.tsv"

exit 0
