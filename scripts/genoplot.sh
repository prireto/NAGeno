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
MUT_LIST="${10}"
BED="${11}"

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

# run for post-filtering data
Rscript ./scripts/depth_analysis.R "$ANALYSIS_DIR/filtered_bam_sr/depth/" "$EXT" "$MOD$MAPQ_MOD" "$ANNO" "$MUT_LIST" "$BED" "$OUT_DIR"

echo "Post-filtering per sample read depth plots saved in: $OUT_DIR as post-filtering_depth.png and post-filtering_depth.svg"

#######################


#echo "############# BAM DEPTH PLOTTING ##############"
echo "This might take a while."

mkdir -p "$ANALYSIS_DIR/filtered_bam_sr/depth/plots"

# run for post-filtering data
ARGS=("$ANALYSIS_DIR/filtered_bam_sr/depth/" "$EXT" "$DEPTH_MOD" "$ANNO" "$MUT_LIST" "$BED")
Rscript "./scripts/depth_analysis.R" "${ARGS[@]}"

echo "Post-filtering per sample read depth plots saved in: $DIR/filtered_bam_sr/depth/plots as post-filtering_depth.svg"

###########################


# run for snp and indel individually
Rscript ./scripts/vcf_amplicon_geno_clairs-to.R "$ANALYSIS_DIR/ClairS-TO/" "$MOD_STRP$MAPQ_MOD$CLAIR_MODEL" "$TXFILE" "$OUT_DIR" "snv" "${SAMPLES[@]}"
Rscript ./scripts/vcf_amplicon_geno_clairs-to.R "$ANALYSIS_DIR/ClairS-TO/" "$MOD_STRP$MAPQ_MOD$CLAIR_MODEL" "$TXFILE" "$OUT_DIR" "indel" "${SAMPLES[@]}"


echo "SNV genotyping info saved in: $OUT_DIR$MOD_STRP$MAPQ_MOD${CLAIR_MODEL}_snv_vcf_collection.tsv."
echo "Indel genotyping info saved in: $OUT_DIR$MOD_STRP$MAPQ_MOD${CLAIR_MODEL}_indel_vcf_collection.tsv."

Rscript ./scripts/genotyping.R "$OUT_DIR" "${MOD_STRP}${MAPQ_MOD}${CLAIR_MODEL}_snv_vcf_collection.tsv" "$MOD_STRP"

echo "Genotyping results saved in: $OUT_DIR$MOD_STRP${MAPQ_MOD}-genotyping_results.tsv"

exit 0
