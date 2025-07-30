#!/bin/bash

# Load parameters

DIR="$1"
ANNO="$2"
REF="$3"
BED="$4"
TXFILE="$5"
THREADS="$6"
MIN_Q="$7"
MAX_U="$8"
MAPQ="$9"
ANALYSIS_DIR="${10}"
EXT="${11}"
CLAIR_PATH="${12}"
CLAIR_MODEL="${13}"
SNPEFF_REF="${14}"

#Automatically assessed
SCRIPT_DIR="${15}"

#preprocessing

mapfile -t SAMPLES < <(cut -f1 $ANNO)
mapfile -t BARCODES < <(cut -f2 $ANNO)

QUANT=$((100 - MAX_U))
MOD="_q${QUANT}_Q${MIN_Q}"
MOD_STRP="q${QUANT}_Q${MIN_Q}"
MAPQ_MOD="_MAPQ${MAPQ}"

DATE=$(date +'%Y-%m-%d')

############################

echo "###############################################################"
echo "Run genotyping on data in $DIR, save in $ANALYSIS_DIR. Use barcodes "${BARCODES[@]}" and corresponding samples "${SAMPLES[@]}". Apply MAPQ filter of $MAPQ and fastq filter of $MOD. As ref use $REF, only run var calling in regions specified in $BED. Var calling model is $CLAIR_MODEL."
echo "###############################################################"

#############################


####################
### filter fastq ###
####################

echo "############# FASTQ FILTERING: $MOD ##############"

echo "Min_Q: $MIN_Q"
echo "MAX_U: $MAX_U"
echo "QUANT: $QUANT"
echo "MOD: $MOD"

mkdir -p "$ANALYSIS_DIR/filtered_fastq"
mkdir -p "$ANALYSIS_DIR/qc/fastplong"

for BC in "${BARCODES[@]}"; do
	echo "${BC}"
	echo "-----"
	FILE="$EXT$BC"
	echo "${FILE}"
	
	# filter fastq
	fastplong --in "$DIR/$FILE.fastq.gz" \
	--qualified_quality_phred $MIN_Q \
	--unqualified_percent_limit $MAX_U \
	--thread $THREADS \
	--out "$ANALYSIS_DIR/filtered_fastq/${FILE}$MOD.fastq.gz" \
	--json "$ANALYSIS_DIR/qc/fastplong/${FILE}$MOD.fastplong.json" \
	--html "$ANALYSIS_DIR/qc/fastplong/${FILE}$MOD.fastplong.html" \
	--disable_adapter_trimming

	echo "${FILE}$MOD.fastq.gz has been created."
done

########################
### align w minimap2 ###
########################

echo "############# ALIGNMENT ##############"


mkdir -p "$ANALYSIS_DIR/filtered_bam_sr"

touch "$ANALYSIS_DIR/filtered_bam_sr/alignment.log"
BAM_LOG="$ANALYSIS_DIR/filtered_bam_sr/alignment.log"

for BC in "${BARCODES[@]}";do

	FILE="$EXT$BC$MOD"

	
	# align fastq
	#use short-read settings (sr) and increase some subsettings like k-mer and minimizer window size
	#generally: increase both for longer reads, use -x map-ont for reads >1kb
	# filter bam file 
	minimap2 -ax sr -k19 -w10 -t "$THREADS" "$REF" "$ANALYSIS_DIR/filtered_fastq/$FILE.fastq.gz" | samtools sort -@ "$THREADS" | samtools view -hbS -q "$MAPQ" -@ "$THREADS" > "$ANALYSIS_DIR/filtered_bam_sr/$FILE$MAPQ_MOD.sorted.bam"

	#index bam files
	samtools index -b "$ANALYSIS_DIR/filtered_bam_sr/$FILE$MAPQ_MOD.sorted.bam" -@ "$THREADS"

	# save alignment stats in log file
	echo "$DATE" >> "$BAM_LOG"
	echo "$FILE" >> "$BAM_LOG"
	samtools flagstat -O tsv -@ "$THREADS" "$ANALYSIS_DIR/filtered_bam_sr/$FILE$MAPQ_MOD.sorted.bam" >> "$BAM_LOG"
	echo "$FILE.fastq has been aligned to $REF and filtered for MAPQ $MAPQ. Stats can be found in $BAM_LOG."
done

################################################
### BAM DEPTH CALC + AMPLICON RANGE PLOTTING ###
################################################

echo "############# BAM DEPTH CALCULATION ##############"

bash "$SCRIPT_DIR/scripts/samtoolsDepth.sh" "$MOD$MAPQ_MOD" "$ANNO" "$BED" "$ANALYSIS_DIR/filtered_bam_sr" "$EXT" "${BARCODES[@]}" "$THREADS"


###################
### VAR CALLING ###
###################

echo "############# VAR CALLING ##############"

for BC in "${BARCODES[@]}";do
	# get rerspective sample name from $ANNO
	SAMPLE=$(awk -v bc="$BC" '$2 == bc {print $1}' "$ANNO")

	FILE="$EXT$BC$MOD$MAPQ_MOD"

	# check if sample name was assigned
	if [[ -n "$SAMPLE" ]]; then
		echo "$BC corresponds to sample $SAMPLE."

		mkdir -p "$ANALYSIS_DIR/ClairS-TO/$SAMPLE$MOD${MAPQ_MOD}_$CLAIR_MODEL"

		# define respective bam file
		BAM="$ANALYSIS_DIR/filtered_bam_sr/$FILE.sorted.bam"

		# run somatic var calling w ClairS-TO
		${CLAIR_PATH}  --tumor_bam_fn "$BAM" --ref_fn "$REF" --threads "$THREADS" --platform "$CLAIR_MODEL" --output_dir "$ANALYSIS_DIR/ClairS-TO/$SAMPLE$MOD${MAPQ_MOD}_$CLAIR_MODEL" --sample_name "$SAMPLE$MOD${MAPQ_MOD}_$CLAIR_MODEL" --snv_min_af 0.05 --indel_min_af 0.1 --min_coverage 4 --qual 12 --python python --samtools samtools --parallel parallel --pypy $(which pypy) --longphase $(which longphase) --whatshap whatshap --bed_fn $BED


		echo "Small somatic variant calling for $FILE is complete using the $CLAIR_MODEL model."
	else
		echo "No corresponding sample name could be found in $ANNO. Update the file and run again."
	fi
done

######################
### vcf annotation ###
######################

echo "############# VCF ANNOTATION ##############"

mkdir -p "$ANALYSIS_DIR/SnpEff"

for SAMPLE in "${SAMPLES[@]}"; do
	# get rerspective sample name from $ANNO
	FILE="$SAMPLE$MOD${MAPQ_MOD}_$CLAIR_MODEL"
	
	# run SnpEff to annotate vcf files
	# allocates more memory and points to the conda/mamba installed .jar file
	java -Xmx8g -jar $(dirname $(which snpEff))/../share/snpeff-*/snpEff.jar -verbose -cancer -lof -stats "$ANALYSIS_DIR/SnpEff/snpEff_summary.html" "$SNPEFF_REF" "$ANALYSIS_DIR/ClairS-TO/$FILE/snv.vcf.gz" > "$ANALYSIS_DIR/ClairS-TO/$FILE/snv.anno.vcf"
	#snpEff -verbose -cancer -lof -stats "$ANALYSIS_DIR/SnpEff/snpEff_summary.html" "$SNPEFF_REF" "$ANALYSIS_DIR/ClairS-TO/$FILE/snv.vcf.gz" > "$ANALYSIS_DIR/ClairS-TO/$FILE/snv.anno.vcf"
	# --canon might be interesting to only include canoncial tx
	# --interval <file> might be interesting to use only regions specified in bed file

	echo "SNV VCF file has been annotated for $FILE using SnpEff and $SNPEFF_REF."
done

######################
### vcf extraction ###
######################

echo "############# VCF EXTRACTION ##############"

# create _temp vcf files w/o headers for all annotated snv files
for SAMPLE in "${SAMPLES[@]}"; do
	# get respective sample name from $ANNO
	FILE="$SAMPLE$MOD${MAPQ_MOD}_$CLAIR_MODEL/snv.anno.vcf"

	# rm all lines until one starts w #CHROM
	sed '/CHROM/,$!d' < "$ANALYSIS_DIR/ClairS-TO/$FILE" > "$ANALYSIS_DIR/ClairS-TO/${FILE}_temp"
	echo "Created ${FILE}_temp."
done

# gunzip all indel files and create _temp vcf files w/o headers
for SAMPLE in "${SAMPLES[@]}"; do
	# get respective sample name from $ANNO
	FILE="$SAMPLE$MOD${MAPQ_MOD}_$CLAIR_MODEL/indel.vcf"

	gunzip "$ANALYSIS_DIR/ClairS-TO/$FILE.gz"

	# rm all lines until one starts w #CHROM
	sed '/CHROM/,$!d' < "$ANALYSIS_DIR/ClairS-TO/$FILE" > "$ANALYSIS_DIR/ClairS-TO/${FILE}_temp"
	echo "Created ${FILE}_temp."
done
