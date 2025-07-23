#!/bin/bash

#save samtools depth of all bam files incl. zero coverage positions
#set thresholds for read quality (-q) or mapping quality (-Q)

#set modifier to subset bam files by (e.g._q90_Q30)
MOD="$1"

BC_ANNO="$2"

BED="$3"

DIR="$4"

OUT_DIR="$DIR/depth"
mkdir -p "$OUT_DIR"

EXT="$5"

THREADS="$6"

mapfile -t BARCODES < <(cut -f2 $BC_ANNO)


######################################

# alt: parallelize for each BAM file up to max thread number
# make "$BED" also accessible for subprocesses
export BED
export OUT_DIR

ls ${DIR}/*${MOD}.sorted.bam | xargs -P "$THREADS" -I {} bash -c '
    BAM="{}"


    echo "Get depth data for $BAM."
    echo -e "contig\tposition\tdepth" > "$OUT_DIR/$(basename "$BAM").depth.tsv"
    samtools depth -b $BED -a "$BAM" >> "$OUT_DIR/$(basename "$BAM").depth.tsv"
    echo "Depth stats saved for $BAM."
'

echo "#############################################################"

######################################

# add actual sample name based on BC

# create an associative array to store barcode to sample name mappings
declare -A barcode_map

echo "Read sample sheet."

# read the barcodes.tsv file line by line
while IFS=$'\t' read -r SAMPLE BC; do
    barcode_map["$BC"]="$SAMPLE"
done < "$BC_ANNO"

echo "Assign sample names to barcodes."

# process each file passed as an argument to the script
for x in "${BARCODES[@]}"; do

	DEPTH_FILE="$OUT_DIR/$EXT$x${MOD}.sorted.bam.depth.tsv"
	echo "$DEPTH_FILE"

    # loop through the barcode mappings and replace the barcodes with sample names
    for BC in "${!barcode_map[@]}"; do

    	if [[ "$BC" == "$x" ]]; then
    		SAMPLE_NAME="${barcode_map[$BC]}"
	        echo "$BC"
	        
	        # use sed to replace the barcode with the sample name in the file
	        mv "$DEPTH_FILE" "$(echo "$DEPTH_FILE" | sed "s/barcode$BC/barcode${BC}_$SAMPLE_NAME/")"
    	fi
    done

    echo "Annotated file: $DEPTH_FILE"
done

echo "All ${MOD}.sorted.bam.depth files processed successfully."

