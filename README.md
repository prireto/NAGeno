
<h3 align="center"> NAGeno - Nanopore Amplicon GENOtyping</h3>

SNV and indel genotyping by Nanopore Amplicon Sequencing and subsequenct variant calling

This pipeline starts with basecalled Nanopore Amplicon sequences and returns 2 overview genotyping tables (SNV and indel), a plot and more elaborate underlying files.
It works for multiplexed samples, as long as each barcode has only been used once.

## Table of contents
* [Introduction](#introduction)
* [Workflow](#workflow)
* [Installation](#installation)
* [Usage](#usage)
* [Tutorial](#tutorial)
* [Citation and Contribution](#citation-and-contribution)
* [License](#license)

## Introduction

What is the issue?

## Workflow
-> Replace with visualisation! -> bring in a picture SHORT description
1) Filter fastq files based on quantile Q scores using **fastplong**: Set the min Q score that a min percentage of the read's bases needs to have (also creates html report) => saves new filtered fastq files
2) Run **multiqc** on fastq files (doesn't seem to work well on the filtered fastqs) => saves multiqc data
3) Align to ref using **minimap2** and sepcial mid-length read settings - saves alignment stats in alignment.log file => saves bam and bai files, and alignment.log
4) Optionally filter bams based on MAPQ score using **samtools** => saves additional bam and bai files
5) calc and plot filtered read depths over amplicons => saves bam.depth.tsv files and depth plots (svg) (for this a bedfile incl. gene names in the 4th col is needed; if several amplicons are in the same gene, make different gene names e.g. gene1_1, gene1_2)
6) SNV and indel variant calling using **clairs-to** (detects germline and somatic mutations even without corresponding normal samples based on SNP databases and specific models) (uses pre-set model, genomic reference, bed file, min quality score and snv/indel min allele frequencies - would be nice if these could be changed as options) => saves one dir per sample and containing several stats and a vcf.gz file for SNVs and indels, respectively
7) vcf file annotation using **SnpEff** (annotates based on genome ref, cancer tag is pre-set here, would be nice to be optional) => returns anno.vcf file
8) Concatenate interesting parts of SNV vcfs into one vcf collection for all samples using **self-written R script** (filters columns of interest but also only the annotation that actually corresponds to one of the panel genes) (needs a tsv file containing a list of genes to filter for in col 1 and the respective transcript id in col 2, should have a col header, is called tx.tsv here) => saves vcf_collection.tsv file
9) Generate overview plots for all identified variants in defined amplicons and one table with only important columns extracted from vcf_collection.tsv using a **self-written R script** => saves a plot and a table

## Installation

Clone this repository

```bash
git clone https://github.com/prinzregententorte/NanoporeAmpliconGenotyping
```

Two `.yaml` files are included into the repository at `envs/scripts`. For the full functionality (i.e. analysis and plotting), both of them need to be created via

```bash
mamba env create -f envs/nageno.yml
mamba env create -f envs/nageno_plot.yml
```

> [!NOTE]
> Alternatively, `conda` can also be used for the creation of the environments, though it will be much slower than using `mamba` or `micromamba`. Since `nageno` ultimately uses the environments for some parts of the analysis, make sure that the command `conda activate nageno` works. Alternatively, make sure that you are setting `--manager` to your respective dependency manager (e.g. `micromamba`) while using the pipeline. 

Further, the somatic variant caller, **ClairS-TO**, and its models need to be installed manually, as explained [here](https://github.com/HKU-BAL/ClairS-TO). 

<details>
<summary>Condensed relevant information about the manual installation of ClairS-TO (click to expand)</summary>

```bash
# create and activate an environment named clairs-to
# install pypy and packages in the environment
# for micromamba
micromamba create -n clairs-to -c bioconda -c pytorch -c conda-forge pytorch tqdm clair3 bcftools einops scipy scikit-learn python=3.9.0 -y
micromamba activate clairs-to

## for anaconda 
#conda create -n clairs-to -c bioconda -c pytorch -c conda-forge pytorch tqdm clair3 bcftools einops python=3.9.0 -y
#source activate clairs-to

git clone https://github.com/HKU-BAL/ClairS-TO.git
cd ClairS-TO

# make sure in clairs-to environment
# download pre-trained models and other resources
echo ${CONDA_PREFIX}
mkdir -p ${CONDA_PREFIX}/bin/clairs-to_models
mkdir -p ${CONDA_PREFIX}/bin/clairs-to_databases
mkdir -p ${CONDA_PREFIX}/bin/clairs-to_cna_data
wget http://www.bio8.cs.hku.hk/clairs-to/models/clairs-to_models.tar.gz
wget http://www.bio8.cs.hku.hk/clairs-to/databases/clairs-to_databases.tar.gz
wget http://www.bio8.cs.hku.hk/clairs-to/cna_data/reference_files.tar.gz
tar -zxvf clairs-to_models.tar.gz -C ${CONDA_PREFIX}/bin/clairs-to_models/
tar -zxvf clairs-to_databases.tar.gz -C ${CONDA_PREFIX}/bin/clairs-to_databases/
tar -zxvf reference_files.tar.gz -C ${CONDA_PREFIX}/bin/clairs-to_cna_data/

./run_clairs_to --help
```

</details>


> [!WARNING]
> `clairs-to` searches for the models at `echo ${CONDA_PREFIX}/bin`. This unfortunately can not be changed easily and thus you need to make sure that `clairs-to_models`, `clairs-to_databases`, and `clairs-to_cna_data` exist in the bin-fodler of the `nageno` environment. You can prevent this extra step by activating the `nageno` environment first and then proceed with the manual `clairs-to` installation. 


## Usage

Generally, `nageno` can be used with two subcommands, `analysis` and `plot`. 

```bash
    Usage: nageno [SUBCOMMAND] [OPTIONS]

Subcommands:

  analysis             Runs genotype analysis. Use --help for mandatory and optional inputs.
  plot                 Runs post-analysis summary and plotting functions.

Use 'nageno [SUBCOMMAND] --help' for more information on a subcommand.

Typical execution order:
  1. Run the analysis subcommand with your parameters:
     nageno analysis [YOUR OPTIONS]

  2. After completion, run the plot subcommand with the same settings:
     nageno plot [YOUR OPTIONS]


```

#### Analysis

```bash
Usage: nageno analysis --dir DIR --anno ANNO --ref REF --bed BED --txfile TXFILE [OPTIONS]

Mandatory arguments:
  --dir                DIR                  Directory containing fastq files
  --anno               ANNO                 Sample sheet file
  --ref                REF                  Reference genome file
  --bed                BED                  BED file for reference
  --txfile             TXFILE               File for visualization

Optional arguments:
  --manager            MANAGER              Package manager used to activate environments (default: conda)
  --threads            THREADS              Number of cores to use (default: 1)
  --min-q              MIN_Q                Minimum base quality (default: 20)
  --max-u              MAX_U                Percentage of bases allowed below MIN_Q (default: 5)
  --mapq               MAPQ                 Minimum mapping quality (default: 0)
  --analysis-dir       DIR                  Directory for output (default: ./analysis)
  --ext                EXT                  Sample name extension (default: SQK-RBK114-24_barcode)
  --clairs-to-path     CLAIR_PATH           Absolute path to 'run_clairs_to'. (default: run_clairs_to)
  --clairs-to-model    CLAIR_MODEL          Clairs-to model (default: ont_r10_dorado_sup_5khz)
  --snpeff-ref         SNPEFF_REF           SNPeff reference genome (default: GRCh38.p14)

```

#### Plot

```bash
Usage: nageno plot --dir DIR --anno ANNO --txfile TXFILE [OPTIONS]

!!! Attention !!!

Make sure you are using the same options as for the analysis.
The generated output files will otherwise not be recognized properly.

Mandatory arguments:
  --dir                DIR                  Directory containing fastq files
  --anno               ANNO                 Sample sheet file
  --txfile             TXFILE               File for visualization
  --bed                BED                  BED file for reference

Optional arguments:
  --manager            MANAGER              Package manager used to activate environments (default: conda)
  --min-q              MIN_Q                Minimum base quality (default: 20)
  --max-u              MAX_U                Percentage of bases allowed below MIN_Q (default: 5)
  --mapq               MAPQ                 Minimum mapping quality (default: 0)
  --analysis-dir       DIR                  Directory for output (default: ./analysis)
  --clairs-to-model    CLAIR_MODEL          Clairs-to model (default: ont_r10_dorado_sup_5khz)
  --mut-list           MUT_LIST             List of mutations for highlighting in depth plots (default: NA)

```

## Tutorial

Using the exemplary test data in `tutorial`, the correct setup can be confirmed and exemplary output can be generated:

```bash
nageno analysis \
  --dir tutorial/data/fastq \
  --anno tutorial/Src/barcode_assignment.tsv \
  --ref /path/to/ref/genome/hg38.fa \
  --bed tutorial/Src/geno_panel_v4.1.bed \
  --txfile ../tutorial/Src/tx.tsv \
  --analysis-dir nageno_tutorial/analysis \ 
  --threads 20 
```

<details>
<summary>Potential installation errors:</summary>

- `[ERROR] file .../envs/nageno/bin/clairs-to_models/ont_r10_dorado_sup_5khz/pileup_affirmative.pkl not found`: Make sure that `clairs-to_models`, `clairs-to_databases`, and `clairs-to_cna_data` exist in the bin-fodler of the `nageno` environment.
- ...

</details>


The `nageno plot` subfunction reulsts in the creation of various different visualisations for the `nageno analysis` output. This is supposed to be used as a quick and comprehensive overview about the genotypes of your samples. 

```bash
nageno plot \
\
\
\
...
```
#### Table 1 –

#### Table 2 –

#### Table 3 – 

#### Plot 1 – Depth plot

#### Plot 2 – Allele frequency plot (protein coding SNVs)

#### Plot 3 – Allele frequency plot (all SNVs)


## Citation and Contribution
BioRXive link / doi

## License

The project is licensed under ...

NOTES: Outdated?
- currently needs to be started NanoporeAmpliconGenotyping dir to work - otherwise it doesn't find the scripts (ERROR: /tmp/tmpyzxjvtjl: line 3: ./scripts/genoSuperscript.sh: No such file or directory)
- maybe as default out-dir create a new dir calles analysis or output in the fastq input location? (mkdir -p $dir/analysis)
- maybe implement a stopping mechanism after an error? or only if 0 files are geenrated in a step? is that too complicated?
- add option for setting absolute clair path
- note in README that calirs-to needs to be installed prior to nagger and provide github link

ERRORS:
- ./scripts/genoSuperscript.sh: line 207: run_clairs_to: command not found
- how about introducing a starting step like --start-at var-calling and internal options to start at any step in the pipeline?

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

