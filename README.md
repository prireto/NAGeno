# NanoporeAmpliconGenotyping
SNV and indel genotyping by Nanopore Amplicon Sequencing and subsequenct variant calling

This pipeline starts with basecalled Nanopore Amplicon sequences and returns 2 overview genotyping table, a plot and more elaborate underlying files.
It works for multiplexed samples, as long as each barcode has only been used once.

**Example cmd:**
bash "/home/vera/gueseer/Scripts/genoSuperscript.sh" 30 10 50 "/home/vera/gueseer/Projects/cancerna/genotyping/np_amplicon_geno/data/amplicon_geno2" "/home/vera/gueseer/Projects/cancerna/genotyping/np_amplicon_geno/analysis/geno2_sr/test" "/home/vera/gueseer/Projects/cancerna/genotyping/np_amplicon_geno/barcode_assignment_geno2.tsv" 01 24

parseable input:
#1: min base Q
#2: percentage that is allowed to miss the minQ (100-x)
#3: MAPQ filter
#4: DIR (data)
#5: ANALYSIS_DIR (output)
#6: ANNO (sample sheet)
#7: BC_START
#8: BC_END



**Workflow:**
1) Filter fastq files based on quantile Q scores using **fastplong**: Set the min Q score that a min percentage of the read's bases needs to have (also creates html repor) => saves new fastq files
2) Run **multiqc** on fastq files (doesn't seem to work well on the filtered fastqs) => saves multiqc data
3) Align to ref using **minimap2** and sepcial mid-length read settings - saves alignment stats in alignment.log file => saves bam and bai files, and alignment.log
4) Optionally filter bams based on MAPQ score using **samtools** => saves additional bam and bai files
( 5) calc and plot depths over amplicons - still needs to be included => saves bam.depth.tsv files and depth plots (svg) )
6) SNV and indel variant calling using **clairs-to** (detects germline and somatic mutations even without corresponding normal samples based on SNP databases and specific models) (uses pre-set model, genomic reference, bed file, min quality score and snv/indel min allele frequencies - would be nice if these could be changed as options) => saves one dir per sample and containing several stats and a vcf.gz file for SNVs and indels, respectively
7) vcf file annotation using **SnpEff** (annotates based on genome ref, cancer tag is pre-set here, would be nice to be optional) => returns anno.vcf file
8) Concatenate interesting parts of vcf into one vcf collection for all samples using **self-written R script** (filters columns of interest but also only the annotation that actually corresponds to one of the panel genes) (needs a tsv file containing a list of genes to filter for in col 1 and the respective transcript id in col 2, should have a col header, is called tx.tsv here) => saves vcf_collection.tsv file
9) Generate overview plots for all identified variants in defined amplicons and one table with only important columns as well as a subset of that one including only protein-coding variants using **self-written R script** (but also synonymous) => saves a plot and 2 tables


###############################################################################################


TODO for pipeline
- in general is it possible to make things as options w default settings?
    - genome ref (for alignment, vcf annotation etc.)
    - kmer and window size for minimap
    - max number of threads to use
    - bed file
    - basically important options for clairs-to
    - nanopore kit
- path to tx.tsv needs to be parsed for vcf extraction - no reasonable default here




Notes
- Pipeline needs to run in clirs-to env
- SnpEff needs newer Java version than we have installed
- SnpEff is installed in /home/vera/gueseer/App/tools/snpEff/


# Background

Genotyping is pain and manual labor. The combination of Oxford Nanopore Technologies (ONT) Amplicon Sequencing and a structured Data Analysis Pipeline can achieve results just as good with less error-prone manual work. 

# Installation and Setup

Nagger can be downloaded from github via .... 

`git clone prinzregententorte/...`

Two .yaml files are included into the repository at envs/scripts. For the full functionality (analysis and plotting), both of them need to be created via

`mamba env create -f envs/nagger.yml`

`mamba env create -f envs/nagger_plotting.yml`

!!! 
    Alternatively, `conda` can also be used for the creation of the environments, though it will be much slower than using `mamba`.

# Usage

Generally, `nagger` can be used with two subcommands, `analysis` and `plotting`. 

```
    ./main.sh 

    Usage: ./main.sh [SUBCOMMAND] [OPTIONS]

    Subcommands:

    analysis             Runs genotype analysis. Use --help for mandatory and optional inputs.
    plot                 Runs post-analysis summary and plotting functions.

    Use './main.sh [SUBCOMMAND] --help' for more information on a subcommand.

    Typical execution order:
    1. Run the analysis subcommand with your parameters:
        ./main.sh analysis [YOUR OPTIONS]

    2. After completion, run the plot subcommand with the same settings:
        ./main.sh plot [YOUR OPTIONS]

```

#### Analysis

```
 ./main.sh analysis --help

Usage: ./main.sh analysis --dir DIR --anno ANNO --ref REF --bed BED --txfile TXFILE [OPTIONS]

Mandatory arguments:
  --dir                DIR                  Directory containing fastq files
  --anno               ANNO                 Sample sheet file
  --ref                REF                  Reference genome file
  --bed                BED                  BED file for reference
  --txfile             TXFILE               File for visualization

Optional arguments:
  --threads            THREADS              Number of cores to use (default: 1)
  --min-q              MIN_Q                Minimum base quality (default: 20)
  --max-u              MAX_U                Percentage of bases allowed below MIN_Q (default: 5)
  --mapq               MAPQ                 Minimum mapping quality (default: 0)
  --analysis-dir       DIR                  Directory for output (default: ./analysis)
  --ext                EXT                  Sample name extension (default: SQK-RBK114-24_barcode)
  --clairs-to-model    CLAIR_MODEL          Clairs-to model (default: _ss)
```

#### Plot

```
./main.sh plot --help

Usage: ./main.sh plot --dir DIR --anno ANNO --txfile TXFILE [OPTIONS]

!!! Attention !!!

Make sure you are using the same options as for the analysis.
The generated output files will otherwise not be recognized properly.

Mandatory arguments:
  --dir                DIR                  Directory containing fastq files
  --anno               ANNO                 Sample sheet file
  --txfile             TXFILE               File for visualization

Optional arguments:
  --min-q              MIN_Q                Minimum base quality (default: 20)
  --max-u              MAX_U                Percentage of bases allowed below MIN_Q (default: 5)
  --mapq               MAPQ                 Minimum mapping quality (default: 0)
  --analysis-dir       DIR                  Directory for output (default: ./analysis)
  --clairs-to-model    CLAIR_MODEL          Clairs-to model (default: _ss)

```

# Visualisation

The `nagger plot` subfunction reulsts in the creation of various different visualisations for the `nagger analysis` output. This is supposed to be used as a quick and comprehensive overview about the genotypes of your samples. 

#### Plot 1

#### Plot 2

#### Plot 3

# License

The project is licensed under ...


NOTES:
- currently needs to be started NanoporeAmpliconGenotyping dir to work - otherwise it doesn't find the scripts (ERROR: /tmp/tmpyzxjvtjl: line 3: ./scripts/genoSuperscript.sh: No such file or directory)
- maybe as default out-dir create a new dir calles analysis or output in the fastq input location? (mkdir -p $dir/analysis)
- maybe implement a stopping mechanism after an error? or only if 0 files are geenrated in a step? is that too complicated?
- add option for setting absolute clair path
- note in README that calirs-to needs to be installed prior to nagger and provide github link

ERRORS:
- ./scripts/genoSuperscript.sh: line 207: run_clairs_to: command not found
- how about introducing a starting step like --start-at var-calling and internal options to start at any step in the pipeline?

