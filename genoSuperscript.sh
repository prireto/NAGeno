#!/bin/bash

# amplicon geno bash superscript
# needs to be in mamba env clairs-to to work!!

################

# HOW TO USE THIS:
# change dir and analysis dir to the dir containing the fastq files of interest and just run w default settings (MAPQ filter=0, q95Q20)

################

# parseable input:
#1: min base Q
#2: percentage that is allowed to miss the minQ (100-x)
#3: MAPQ filter
#4: DIR (data)
#5: ANALYSIS_DIR (output)
#6: ANNO (sample sheet)
#7: BC_START
#8: BC_END


################
### SET VARs ###
################

# fastq dir
#DIR="/home/vera/gueseer/Projects/cancerna/genotyping/np_amplicon_geno/data/amplicon_geno2"
DIR=$4

#ANALYSIS_DIR="/home/vera/gueseer/Projects/cancerna/genotyping/np_amplicon_geno/analysis/geno2_sr"
ANALYSIS_DIR=$5


# samplesheet
#ANNO="/home/vera/gueseer/Projects/cancerna/genotyping/np_amplicon_geno/barcode_assignment.tsv"
ANNO=$6

# get sample names from sample sheet (first col, no header) adn store in array called $SAMPLES
mapfile -t SAMPLES < <(cut -f1 $ANNO)
# same for barcodes in $BARCODES - could replace {01..24} by "${SAMPLES[@]}"
mapfile -t BARCODES < <(cut -f2 $ANNO)

# Bed file for geno seq ref
BED="/home/vera/gueseer/Src/geno_panel_v2.1.bed"

# Check if input arguments are provided for BC range, otherwise use defaults
#if [[ -z "$7" && -z "$8" ]]; then
#	# use pre-sets
#	BC_START="$DEFAULT_BC_START"
#	BC_END="$DEFAULT_BC_END"
#	echo "No values parsed for barcode range."
#else
#	# use parsed values
#	BC_START=$7
#	BC_END=$8
#	echo "Use parsed values for barcode range."
#fi

#create barcode array - if we do it like this this could also just be a parseable array in case it's not a consecutive list of BCs
# BARCODES=$(eval echo {$BC_START..$BC_END})


# sample name extension
EXT="SQK-RBK114-24_barcode"

# Modifier based on filtering - e.g. _q90_Q20
#MOD="_q90_Q30"

#if values are parsed use these, otherwise use default
#set quantile and Q threshold for fastplong fastq filtering = default values
# also set MAPQ filtering threshold
# Default values
DEFAULT_MIN_Q=20
DEFAULT_MAX_U=5
DEFAULT_MAPQ=0

# Check if input arguments are provided, otherwise use defaults
if [[ -z "$1" && -z "$2" && -z "$3" ]]; then
	# use pre-sets
	MIN_Q="$DEFAULT_MIN_Q"
	MAX_U="$DEFAULT_MAX_U"
	MAPQ="$DEFAULT_MAPQ"
	echo "No values parsed for fastq and bam filtering."
else
	# use parsed values
	MIN_Q=$1
	MAX_U=$2
	MAPQ=$3
	echo "Use parsed values for fastq and bam filtering."
fi

QUANT=$(( 100-MAX_U ))
MOD="_q${QUANT}_Q${MIN_Q}"
MOD_STRP="q${QUANT}_Q${MIN_Q}"

# ref fa file - here just geno panel ref (must be samtools indexed)
REF_GENO="/home/vera/gueseer/Src/GRCh38.p14.genome.geno-panel.fa"
REF_WHOLE="/home/vera/gueseer/Src/GRCh38.p14.genome.fa"
REF_CHR="/home/vera/gueseer/Src/GRCh38.p14.genome.chr.fa"

#set ref
REF="$REF_WHOLE"

# Bed file for geno seq ref
BED="/home/vera/gueseer/Src/geno_panel_v2.1.bed"
# seq_region bed (based on samtools coverage decision, not primer_design)
#BED="/home/vera/gueseer/Src/geno_panel_seq_region_v1.2.bed"

# set model used for var calling
# ssrs is a model trained initially with synthetic samples and then real samples augmented (e.g., ont_r10_dorado_sup_5khz_ssrs), ss is a model trained from synthetic samples (e.g., ont_r10_dorado_sup_5khz_ss). The ssrs model provides better performance and fits most usage scenarios. ss model can be used when missing a cancer-type in model training is a concern. In v0.3.0, four real cancer cell-line datasets (HCC1937, HCC1954, H1437, and H2009) covering two cancer types (breast cancer, lung cancer) published by Park et al. were used for ssrs model training.

CLAIR_MODEL="_ss"


DATE=$(date +'%Y-%m-%d')

##
# also check barcodes and dir names in following code!

############################

echo "###############################################################"
echo "Run genotyping on data in $DIR, save in $ANALYSIS_DIR. Use barcodes "${BARCODES[@]}" and corresponding samples "${SAMPLES[@]}". Apply MAPQ filter of $MAPQ and fastq filter of $MOD. As ref use $REF, only run var calling in regions specified in $BED. Var calling model is $CLAIR_MODEL."
echo "###############################################################"

#############################

# optionally include basecalling here?


####################
### filter fastq ###
####################

echo "############# FASTQ FILTERING: $MOD ##############"

echo "Min_Q: $MIN_Q"
echo "MAX_U: $MAX_U"
echo "QUANT: $QUANT"
echo "MOD: $MOD"

# q90 Q20 (90% of bases in a read need to have base q of at least 20)

# create dir if it doesn't exist yet
	mkdir -p "$DIR/filtered_fastq"

for BC in "${BARCODES[@]}"; do

	FILE="$EXT$BC"

	# if out file already exists skip
	if [[ -e  "$DIR/filtered_fastq/${FILE}$MOD.fastq.gz" ]]; then
		# skip
		echo "${FILE}$MOD.fastq.gz already exists."

	else
		# filter fastq
		fastplong -i "$DIR/fastq/$FILE.fastq.gz" --qualified_quality_phred $MIN_Q --unqualified_percent_limit $MAX_U --thread 16 --out "$DIR/filtered_fastq/${FILE}$MOD.fastq.gz" --json "$DIR/filtered_fastq/fastplong.json" --html "$DIR/filtered_fastq/fastplong.html"

		echo "${FILE}$MOD.fastq.gz has been created."
	fi
done

# include multiqc? pre and post filtering?

###############
### multiqc ###
###############

echo "🎄"

multiqc --force --dirs "$DIR/filtered_fastq" --outdir "$DIR/filtered_fastq" --verbose

echo "multiqc performed on all filtered samples."


########################
### align w minimap2 ###
########################

# based on minimap2-midLengthBC.sh

echo "############# ALIGNMENT ##############"

# create dir if it doesn't exist yet
	mkdir -p "$DIR/filtered_bam_sr"

# create alignment log file if it doesn't exist yet
touch "$DIR/filtered_bam_sr/alignment.log"
BAM_LOG="$DIR/filtered_bam_sr/alignment.log"

for BC in "${BARCODES[@]}";do

	FILE="$EXT$BC$MOD"

	# skip if out file already exists
	if [[ -e  "$DIR/filtered_bam_sr/$FILE.sorted.bam" ]]; then
		# skip
		echo "$FILE.sorted.bam already exists."

		#echo "$FILE"
		#echo "$DIR/fastq/$FILE.fastq.gz"
		#echo "$DIR/bam_sr/$FILE.sorted.bam"

	else
		# align fastq
		#use short-read settings (sr) and increase some subsettings like k-mer and minimizer window size
		#generally: increase both for longer reads, use -x map-ont for reads >1kb
		#kmer size for breaking down reads for indexing and alignment
		#lower window sides means less comparisons => less sensitivity for quicker results
		minimap2 -ax sr -k19 -w10 -t '60' "$REF" "$DIR/filtered_fastq/$FILE.fastq.gz" | samtools sort -@ 60 | samtools view -hbS -@ 60 > "$DIR/filtered_bam_sr/$FILE.sorted.bam"

		#index bam files
		samtools index -b "$DIR/filtered_bam_sr/$FILE.sorted.bam" -@ 60


		# save alignment stats in log file
		echo "$DATE" >> "$BAM_LOG"
		echo "$FILE" >> "$BAM_LOG"
		samtools flagstat -O tsv -@ 60 "$DIR/filtered_bam_sr/$FILE.sorted.bam" >> "$BAM_LOG"
		echo "$FILE.fastq has been aligned to $REF. Stats can be found in $BAM_LOG."
	fi

	# optionally filter bam file (usually for MAPQ>50)
	MAPQ_MOD="_MAPQ${MAPQ}"
	if [[ "$MAPQ" -ne 0 ]]; then
		#set MAPQ filter var (do filtering or not)
		MAPQ_FILTER=1
	else
		MAPQ_FILTER=0
	fi
	# skip if out file already exists
	if [[ ! -e  "$DIR/filtered_bam_sr/$FILE$MAPQ_MOD.sorted.bam" && "$MAPQ_FILTER" -eq 1 ]]; then
		echo "Filter bams by MAPQ$MAPQ."

		samtools view -b -q 50 -@ 60 "$DIR/filtered_bam_sr/$FILE.sorted.bam" > "$DIR/filtered_bam_sr/$FILE$MAPQ_MOD.sorted.bam"

		#index MAPQ filtered bam files
		samtools index -b "$DIR/filtered_bam_sr/$FILE$MAPQ_MOD.sorted.bam" -@ 60

		echo "$FILE.sorted.bam has been filtered for MAPQ$MAPQ."
	fi
done


# BAM refinement - base q score recalibration?
# duplicate removal? doesn't really make sense here bc all fragments should be amplicons of same region


################################################
### BAM DEPTH CALC + AMPLICON RANGE PLOTTING ###
################################################

# use samtoolsDepth.sh


###################
### VAR CALLING ###
###################

echo "############# VAR CALLING ##############"

# use somatic small var caller clairS TO
# needs to be in mamba env clairs-to to work!!

for BC in "${BARCODES[@]}";do
	# get rerspective sample name from $ANNO
	SAMPLE=$(awk -v bc="$BC" '$2 == bc {print $1}' "$ANNO")

	FILE="$EXT$BC$MOD$MAPQ_MOD"

	# check if sample name was assigned
	if [[ -n "$SAMPLE" ]]; then
		# proceed
		echo "$BC corresponds to sample $SAMPLE."

		# skip if out dir already exists
		if [[ -e  "$ANALYSIS_DIR/ClairS-TO/$SAMPLE$MOD$MAPQ_MOD$CLAIR_MODEL" ]]; then
			# skip
			echo "$ANALYSIS_DIR/ClairS-TO/$SAMPLE$MOD$MAPQ_MOD$CLAIR_MODEL already exists. Apparently small somatic variant calling has already been performed."
		else
			# create out dir
			mkdir -p "$ANALYSIS_DIR/ClairS-TO/$SAMPLE$MOD$MAPQ_MOD$CLAIR_MODEL"

			# define respective bam file
			BAM="$DIR/filtered_bam_sr/$FILE.sorted.bam"

			# run somatic var calling w ClairS-TO
			/home/vera/gueseer/App/tools/ClairS-TO/run_clairs_to  --tumor_bam_fn "$BAM" --ref_fn "$REF" --threads 60 --platform "ont_r10_dorado_sup_5khz$CLAIR_MODEL" --output_dir "$ANALYSIS_DIR/ClairS-TO/$SAMPLE$MOD$MAPQ_MOD$CLAIR_MODEL" --sample_name "$SAMPLE$MOD$MAPQ_MOD$CLAIR_MODEL" --snv_min_af 0.05 --indel_min_af 0.1 --min_coverage 4 --qual 12 --python "/home/vera/gueseer/App/mamba/envs/clairs-to/bin/python" --pypy "/home/vera/gueseer/App/mamba/envs/clairs-to/bin/pypy3" --samtools "/home/vera/gueseer/App/mamba/envs/clairs-to/bin/samtools" --parallel "/home/vera/gueseer/App/mamba/envs/clairs-to/bin/parallel" --longphase "/home/vera/gueseer/App/mamba/envs/clairs-to/bin/longphase" --whatshap "/home/vera/gueseer/App/mamba/envs/clairs-to/bin/whatshap" --bed_fn $BED


			echo "Small somatic variant calling for $FILE is complete using the $CLAIR_MODEL model."
		fi
	else
		echo "No corresponding sample name could be found in $ANNO. Update the file and run again."
	fi
done
# also try with ont_r10_dorado_sup_5khz_ss
# MODEL="r1041_e82_400bps_sup"
# --genotyping_mode_vcf_fn VCF_FILE for only genotyping specific hotspots
#  --bed_fn "$BED"
# --region chr2:197389784-197435093

# post hoc filtering of vcf file to select only "good" vars?
# vcf processing w bgzip and tabix?

######################
### vcf annotation ###
######################

echo "############# VCF ANNOTATION ##############"


# initially wanted to use vep but there are too many issues - use SnpEff instead
# filter and annotate vcfs using ensemble-vep (see /home/vera/gueseer/Scripts/vep.sh or vepDocker.sh)
# http://www.ensembl.org/info/docs/tools/vep/script/vep_options.html

#ONLY WORKS with more recent java version than globally installed one - use the one from snpeff env

mkdir "$ANALYSIS_DIR/SnpEff"

for SAMPLE in "${SAMPLES[@]}"; do
	# get rerspective sample name from $ANNO
	FILE="$SAMPLE$MOD$MAPQ_MOD$CLAIR_MODEL"
	SNPEFF_REF="GRCh38.p14"

	# skip if out file already exists
	if [[ -e  "$ANALYSIS_DIR/ClairS-TO/$FILE/snv.anno.vcf" ]]; then
		# skip
		echo "$ANALYSIS_DIR/ClairS-TO/$FILE/snv.anno.vcf already exists. Apparently vcf annotation has already been performed."
	else
		# run SnpEff to annotate vcf files
		# vcf.gz works as input, output is uncompressed vcf
		/home/vera/gueseer/App/mamba/envs/snpeff/bin/java -Xmx8g -jar /home/vera/gueseer/App/tools/snpEff/snpEff.jar -verbose -cancer -lof -stats "$ANALYSIS_DIR/SnpEff/snpEff_summary.html" "$SNPEFF_REF" "$ANALYSIS_DIR/ClairS-TO/$FILE/snv.vcf.gz" > "$ANALYSIS_DIR/ClairS-TO/$FILE/snv.anno.vcf"

		# --canon might be interesting to only include canoncial tx
		# --interval <file> might be interesting to use only regions specified in bed file

		echo "SNV VCF file has been annotated for $FILE using SnpEff and $SNPEFF_REF."
	fi
done



######################
### vcf extraction ###
######################

echo "############# VCF EXTRACTION ##############"

# create _temp vcf files w/o headers
for SAMPLE in "${SAMPLES[@]}"; do
	# get respective sample name from $ANNO
	FILE="$SAMPLE$MOD$MAPQ_MOD$CLAIR_MODEL/snv.anno.vcf"

	# rm all lines until one starts w #CHROM
	sed '/CHROM/,$!d' < "$ANALYSIS_DIR/ClairS-TO/$FILE" > "$ANALYSIS_DIR/ClairS-TO/${FILE}_temp"
	echo "Created ${FILE}_temp."
done


#run Rscript431 vcf_amplicon_geno.R providing all samples to process all of the vcf files for these samples and create a common table for all samples (Date_vcf_collection.tsv)

#############################################
###### what the R script actually does ######
#############################################

# sep cols by TAB
# mv both last cols before FORMAT
# sep "SAMPLE" col by ':', take col names by separating FORMAT col by ':'
# delete FORMAT col
# sep cols based on ';' - generic names? or INFOx
# rm all cols after FORMAT col
# sep cols based on ','
# rm all post FORMAT cols except the one that matches tx#
# sep cols based on '|'
# more stuff - all explained in detail in R script
##############################################

# expand array into separate args to pass to Rscript ("${SAMPLES[@]}") - otherwise they are not properly embedded into a String array
# Rscript needs (1) dir w a folder per sample containing the vcf files (2) $MOD$MAPQ_MOD$CLAIR_MODEL (3) all samples as String array => all in one collective String array
#create the String array to parse to Rscript
ARGS=("$ANALYSIS_DIR/ClairS-TO/" "$MOD_STRP$MAPQ_MOD$CLAIR_MODEL" "${SAMPLES[@]}")
/home/vera/eberhamn/lib/R-4.3.1/bin/Rscript "/home/vera/gueseer/Pipelines/veralab_geno_vcf/vcf_amplicon_geno_clairs-to.R" "${ARGS[@]}"


# can there be an arg to make this optional?
# remove vcf temp files again
#for SAMPLE in "${SAMPLES[@]}"; do
#	rm "$ANALYSIS_DIR/ClairS-TO/${FILE}_temp"
#	echo "Removed ${FILE}_temp."
#
#	echo "Information relevant for genotyping has been extracted from $FILE."
#	echo "######################"
#done


echo "A table with all the genotyping info has been created: $MOD_STRP$MAPQ_MOD${CLAIR_MODEL}_vcf_collection.tsv."

##########################################################################

# add step that emits table and plot for genotyped samples 

ARGS2=("$ANALYSIS_DIR/ClairS-TO/" "$MOD_STRP$MAPQ_MOD${CLAIR_MODEL}_vcf_collection.tsv" "$MOD_STRP")
/home/vera/eberhamn/lib/R-4.3.1/bin/Rscript "/home/vera/gueseer/Scripts/genotyping.R" "${ARGS2[@]}"

echo "Genotyping results saved in $ANALYSIS_DIR/ClairS-TO/${MOD_STRP}-genotyping_results.tsv, -genotyping_prot_coding_results.tsv and prot coding muts plotted in ${MOD_STRP}-genotyping_results.svg"


#save this file as log in var calling results dir
cp $0 "$ANALYSIS_DIR/ClairS-TO/"
