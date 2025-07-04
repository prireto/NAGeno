
<h3 align="center"> NAGeno - Nanopore Amplicon GENOtyping</h3>

A comprehensive pipeline for SNV and indel genotyping on Nanopore Amplicon Sequencing data.

NAGeno starts with basecalled Nanopore Amplicon sequences and returns two overview genotyping tables (SNV and indel), a SNV genotype overview plot and more elaborate underlying files.
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

**Accurate genotyping made simple.**

Identifying SNVs and indels is essential in molecular biology and clinical diagnostics. Sangersequencing—still the gold-standard for its high accuracy—requires manual inspection to avoid artifacts and catch low-frequency variants. However, it often struggles in GC-rich or highly repetitive regions. NGS provides even higher accuracy with automated analysis, but is typically excessive for small to medium-scale projects and routine lab workflows.

**NAGeno** - **N**anopore **A**mplicon **Geno**typing combines high accuracy even in GC and reptitive regions with Sanger-like simplicity while ensuring scalability making genotyping both robust and effortless.

## Workflow

NAGeno performs SNV and indel genotyping on fastq files of nanopore amplicon sequencing. Amplicons can cover regions of approx. 50 bp - 5 kb. In a fully automated workflow, we generate detailed tables for both SNVs and indels, along with an overview plot of SNVs per sample, by following these steps:

<div align="center">
    <img src="https://github.com/user-attachments/assets/6ef939a9-2c55-47db-b308-c4d900f81268" width="400">
</div>

## Installation

Clone this repository

```bash
git clone https://github.com/prinzregententorte/NanoporeAmpliconGenotyping
```

Two `.yaml` files are included into the repository at `envs/scripts`. For the full functionality (i.e. analysis and plotting), both of them need to be created via

```bash
conda env create -f NanoporeAmpliconGenotyping/envs/nageno.yml
conda env create -f NanoporeAmpliconGenotyping/envs/nageno_plot.yml
conda activate nageno
```

> [!NOTE]
> Alternative to `conda`, `mamba` or `micromamba` can also be used for the creation of the environments which will be much faster. Since `nageno` ultimately uses the environments for some parts of the analysis, make sure that the command `conda activate nageno` works. Alternatively, make sure that you are setting `--manager` to your respective dependency manager (e.g. `micromamba`) while using the pipeline. 

Further, the somatic variant caller, **ClairS-TO**, and its models need to be installed manually, as explained [here](https://github.com/HKU-BAL/ClairS-TO). 

<details>
<summary>Condensed relevant information about the manual installation of ClairS-TO (click to expand)</summary>

```bash
#SRP DELETE?
# create and activate an environment named clairs-to
# install pypy and packages in the environment

# for mamba
#mamba create -n clairs-to -c bioconda -c pytorch -c conda-forge pytorch tqdm clair3 bcftools einops scipy scikit-learn python=3.9.0 -y
#mamba activate clairs-to

# for micromamba
#micromamba create -n clairs-to -c bioconda -c pytorch -c conda-forge pytorch tqdm clair3 bcftools einops scipy scikit-learn python=3.9.0 -y
#micromamba activate clairs-to

## for anaconda 
#conda create -n clairs-to -c bioconda -c pytorch -c conda-forge pytorch tqdm clair3 bcftools einops python=3.9.0 -y
#source activate clairs-to

# in case of a timeout error (Download error (28) Timeout was reached) try modifying timeout settings (works like this only for mamba and conda)
#conda config --set remote_connect_timeout_secs 30
#conda config --set remote_read_timeout_secs 30

#SRP DELETE? END
 

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
> `clairs-to` searches for the models at `echo ${CONDA_PREFIX}/bin`. This unfortunately can not be changed easily and thus you need to make sure that `clairs-to_models`, `clairs-to_databases`, and `clairs-to_cna_data` exist in the bin-folder of the `nageno` environment. You can prevent this extra step by, as described above, activating the `nageno` environment first and then proceed with the manual `clairs-to` installation.

SRP: if you want to access the tool from anywhere not just the dir you installed it in, you can add the path to you bashrc like this:
XXXX
Actually, you need to actually be at the nageno place because all paths are relative to that




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
  --txfile tutorial/Src/tx.tsv \
  --analysis-dir tutorial/analysis \ 
  --threads 20 
```

<details>
<summary>Potential installation errors:</summary>

- `[ERROR] file .../envs/nageno/bin/clairs-to_models/ont_r10_dorado_sup_5khz/pileup_affirmative.pkl not found`: Make sure that `clairs-to_models`, `clairs-to_databases`, and `clairs-to_cna_data` exist in the bin-folder of the `nageno` environment.
- ...

</details>


The `nageno plot` subfunction reulsts in the creation of various different visualisations for the `nageno analysis` output. This is supposed to be used as a quick and comprehensive overview about the genotypes of your samples.

```bash
nageno plot \
  --dir tutorial/data/fastq \
  --anno tutorial/Src/barcode_assignment.tsv \
  --ref /path/to/ref/genome/hg38.fa \
  --bed tutorial/Src/geno_panel_v4.1.bed \
  --txfile ../tutorial/Src/tx.tsv \
  --analysis-dir nageno_tutorial/analysis \ 
  --threads 20 
```

> [!TIP]
> `nageno plot` needs less arguments than `nageno analysis`. Since additional arguments are ignored, the quickest way to use the plotting functionality on your results is by replacing the `analysis` with the `plot` subcommand and re-run.  

#### Table 1 –

#### Table 2 –

#### Table 3 – 

#### Plot 1 – Depth plot

#### Plot 2 – Allele frequency plot (protein coding SNVs)

#### Plot 3 – Allele frequency plot (all SNVs)


## Citation and Contribution
BioRXive link / doi

## License

This project is licensed under the [Apache License 2.0](LICENSE).
[![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg)](LICENSE)

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

