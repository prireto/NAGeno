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

Rscript ./scripts/vcf_amplicon_geno_clairs-to.R "$ANALYSIS_DIR/ClairS-TO/" "$MOD_STRP$MAPQ_MOD$CLAIR_MODEL" "$TXFILE" "$OUT_DIR" "${SAMPLES[@]}"

echo "Genotyping info saved in: $OUT_DIR$MOD_STRP$MAPQ_MOD$CLAIR_MODEL_vcf_collection.tsv"

Rscript ./scripts/genotyping.R "$OUT_DIR" "${MOD_STRP}${MAPQ_MOD}${CLAIR_MODEL}_vcf_collection.tsv" "$MOD_STRP"

echo "Genotyping results saved in: $OUT_DIR/${MOD_STRP}-genotyping_results.tsv"

exit 0