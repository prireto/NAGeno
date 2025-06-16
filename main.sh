#!/bin/bash
# Default values
DEFAULT_THREADS=1
DEFAULT_MIN_Q=20
DEFAULT_MAX_U=5
DEFAULT_MAPQ=0
DEFAULT_ANALYSIS_DIR="./analysis"
DEFAULT_EXT="SQK-RBK114-24_barcode"
DEFAULT_CLAIR_PATH="run_clairs_to"
DEFAULT_CLAIR_MODEL="_ss"
DEFAULT_PLOT_ONLY=0
DEFAULT_MUT_LIST="NA" # I dont know if this is used correctly as NA in R now. Needs to be checked in the R script.

# Help function for main script
display_main_help() {
    echo
    echo "Usage: $0 [SUBCOMMAND] [OPTIONS]"
    echo
    echo "Subcommands:"
    echo
    printf "  %-20s %s\n" "analysis" "Runs genotype analysis. Use --help for mandatory and optional inputs."
    printf "  %-20s %s\n" "plot" "Runs post-analysis summary and plotting functions."
    echo
    echo "Use '$0 [SUBCOMMAND] --help' for more information on a subcommand."
    echo
    echo "Typical execution order:"
    echo "  1. Run the analysis subcommand with your parameters:"
    echo "     $0 analysis [YOUR OPTIONS]"
    echo
    echo "  2. After completion, run the plot subcommand with the same settings:"
    echo "     $0 plot [YOUR OPTIONS]"
    echo
    exit 0
}

# Help function for plot subcommand
display_plot_help() {
    echo
    echo "Usage: $0 plot --dir DIR --anno ANNO --txfile TXFILE [OPTIONS]"
    echo
    echo "!!! Attention !!!"
    echo
    echo -e "Make sure you are using the same options as for the analysis.\nThe generated output files will otherwise not be recognized properly."
    echo
    echo "Mandatory arguments:"
    printf "  %-20s %-20s %s\n" "--dir" "DIR" "Directory containing fastq files"
    printf "  %-20s %-20s %s\n" "--anno" "ANNO" "Sample sheet file"
    printf "  %-20s %-20s %s\n" "--txfile" "TXFILE" "File for visualization"
    printf "  %-20s %-20s %s\n" "--bed" "BED" "BED file for reference"
    echo
    echo "Optional arguments:"
    printf "  %-20s %-20s %s\n" "--min-q" "MIN_Q" "Minimum base quality (default: 20)"
    printf "  %-20s %-20s %s\n" "--max-u" "MAX_U" "Percentage of bases allowed below MIN_Q (default: 5)"
    printf "  %-20s %-20s %s\n" "--mapq" "MAPQ" "Minimum mapping quality (default: 0)"
    printf "  %-20s %-20s %s\n" "--analysis-dir" "DIR" "Directory for output (default: ./analysis)"
    printf "  %-20s %-20s %s\n" "--clairs-to-model" "CLAIR_MODEL" "Clairs-to model (default: _ss)"
    printf "  %-20s %-20s %s\n" "--mut-list" "MUT_LIST" "List of mutations for highlighting in depth plots (default: NA)"
    echo
    exit 0
}

# Help function for analysis subcommand
display_analysis_help() {
    echo
    echo "Usage: $0 analysis --dir DIR --anno ANNO --ref REF --bed BED --txfile TXFILE [OPTIONS]"
    echo
    echo "Mandatory arguments:"
    printf "  %-20s %-20s %s\n" "--dir" "DIR" "Directory containing fastq files"
    printf "  %-20s %-20s %s\n" "--anno" "ANNO" "Sample sheet file"
    printf "  %-20s %-20s %s\n" "--ref" "REF" "Reference genome file"
    printf "  %-20s %-20s %s\n" "--bed" "BED" "BED file for reference"
    printf "  %-20s %-20s %s\n" "--txfile" "TXFILE" "File for visualization"
    echo
    echo "Optional arguments:"
    printf "  %-20s %-20s %s\n" "--threads" "THREADS" "Number of cores to use (default: 1)"
    printf "  %-20s %-20s %s\n" "--min-q" "MIN_Q" "Minimum base quality (default: 20)"
    printf "  %-20s %-20s %s\n" "--max-u" "MAX_U" "Percentage of bases allowed below MIN_Q (default: 5)"
    printf "  %-20s %-20s %s\n" "--mapq" "MAPQ" "Minimum mapping quality (default: 0)"
    printf "  %-20s %-20s %s\n" "--analysis-dir" "DIR" "Directory for output (default: ./analysis)"
    printf "  %-20s %-20s %s\n" "--ext" "EXT" "Sample name extension (default: SQK-RBK114-24_barcode)"
    printf "  %-20s %-20s %s\n" "--clairs-to-path" "CLAIR_PATH" "Absolute path to 'run_clairs_to'. (default: run_clairs_to)"
    printf "  %-20s %-20s %s\n" "--clairs-to-model" "CLAIR_MODEL" "Clairs-to model (default: _ss)"
    echo
    exit 0
}

# Check for subcommands and help options
if [[ $# -eq 0 ]]; then
    display_main_help
fi

SUBCOMMAND="$1"
shift

case "$SUBCOMMAND" in
    plot)
        if [[ "$1" == "--help" ]]; then
            display_plot_help
        fi
        ;;
    analysis)
        if [[ "$1" == "--help" ]]; then
            display_analysis_help
        fi
        ;;
    -h|--help)
        display_main_help
        ;;
    *)
        echo "Unknown subcommand: $SUBCOMMAND"
        display_main_help
        ;;
esac

# Parse arguments for plot and analysis subcommands
while [[ $# -gt 0 ]]; do
    case "$1" in
        --dir) DIR="$2"; shift 2;;
        --anno) ANNO="$2"; shift 2;;
        --ref) REF="$2"; shift 2;;
        --bed) BED="$2"; shift 2;;
        --txfile) TXFILE="$2"; shift 2;;
        --threads) THREADS="$2"; shift 2;;
        --min-q) MIN_Q="$2"; shift 2;;
        --max-u) MAX_U="$2"; shift 2;;
        --mapq) MAPQ="$2"; shift 2;;
        --analysis-dir) ANALYSIS_DIR="$2"; shift 2;;
        --ext) EXT="$2"; shift 2;;
        --clairs-to-path) CLAIR_PATH="$2"; shift 2;;
        --clairs-to-model) CLAIR_MODEL="$2"; shift 2;;
        --plot-only) PLOT_ONLY=1; shift;;
        --mut-list) MUT_LIST="$2"; shift 2;;
        *) echo "Unknown option: $1"; display_main_help;;
    esac
done

# Check mandatory arguments for plot subcommand
if [[ "$SUBCOMMAND" == "plot" && ( -z "$DIR" || -z "$ANNO" || -z "$TXFILE" ) ]]; then
    echo "Error: Missing mandatory argument(s) for plot. Use --help for usage."
    exit 1
fi

# Check mandatory arguments for analysis subcommand
if [[ "$SUBCOMMAND" == "analysis" && ( -z "$DIR" || -z "$ANNO" || -z "$REF" || -z "$BED" || -z "$TXFILE" ) ]]; then
    echo "Error: Missing mandatory argument(s) for analysis. Use --help for usage."
    exit 1
fi

# Assign default values if not provided
THREADS=${THREADS:-$DEFAULT_THREADS}
MIN_Q=${MIN_Q:-$DEFAULT_MIN_Q}
MAX_U=${MAX_U:-$DEFAULT_MAX_U}
MAPQ=${MAPQ:-$DEFAULT_MAPQ}
ANALYSIS_DIR=${ANALYSIS_DIR:-$DEFAULT_ANALYSIS_DIR}
EXT=${EXT:-$DEFAULT_EXT}
CLAIR_PATH=${CLAIR_PATH:-$DEFAULT_CLAIR_PATH}
CLAIR_MODEL=${CLAIR_MODEL:-$DEFAULT_CLAIR_MODEL}
PLOT_ONLY=${PLOT_ONLY:-$DEFAULT_PLOT_ONLY}
MUT_LIST=${MUT_LIST:-$DEFAULT_MUT_LIST}

# Print parsed values
echo "Using the following settings:"
echo "  THREADS=$THREADS"
echo "  MIN_Q=$MIN_Q"
echo "  MAX_U=$MAX_U"
echo "  MAPQ=$MAPQ"
echo "  DIR=$DIR"
echo "  ANALYSIS_DIR=$ANALYSIS_DIR"
echo "  ANNO=$ANNO"
echo "  REF=$REF"
echo "  BED=$BED"
echo "  TXFILE=$TXFILE"
echo "  EXT=$EXT"
echo "  CLAIR_PATH=$CLAIR_PATH"
echo "  CLAIR_MODEL=$CLAIR_MODEL"
echo "  PLOT_ONLY=$PLOT_ONLY"
echo "  MUT_LIST=$MUT_LIST"

# Parameter preprocessing
QUANT=$((100 - MAX_U))
MOD="_q${QUANT}_Q${MIN_Q}"
MOD_STRP="q${QUANT}_Q${MIN_Q}"
MAPQ_MOD="_MAPQ${MAPQ}"
DATE=$(date +'%Y-%m-%d')

# log files
LOG_DIR="$ANALYSIS_DIR/logs"
mkdir -p "$LOG_DIR"  # Ensure log directory exists

TIMESTAMP=$(date +"%Y-%m-%d_%H-%M-%S")

if [[ "$SUBCOMMAND" == "plot" ]]; then
    LOG_FILE="$LOG_DIR/plot_${TIMESTAMP}.log"
    
    echo "Running plot script... Logging to $LOG_FILE"
    
    # Run the script and log output
    conda run --no-capture-output -n nagger_plotting \
        ./scripts/genoplot.sh "$DIR" "$ANNO" "$TXFILE" "$ANALYSIS_DIR" "$CLAIR_MODEL" "$MIN_Q" "$MAX_U" "$MAPQ" "$EXT" "$MUT_LIST" "$BED"\
        | tee "$LOG_FILE"

    echo "Plot script finished! Log saved at $LOG_FILE"

elif [[ "$SUBCOMMAND" == "analysis" ]]; then
    LOG_FILE="$LOG_DIR/analysis_${TIMESTAMP}.log"

    echo "Running analysis script... Logging to $LOG_FILE"

    # Run the script and log output
    conda run --no-capture-output -n nagger \
        ./scripts/genoSuperscript.sh "$DIR" "$ANNO" "$REF" "$BED" "$TXFILE" "$THREADS" "$MIN_Q" "$MAX_U" "$MAPQ" "$ANALYSIS_DIR" "$EXT" "$CLAIR_PATH" "$CLAIR_MODEL" \
        | tee "$LOG_FILE"

    echo "Analysis script completed! Log saved at $LOG_FILE"
fi
