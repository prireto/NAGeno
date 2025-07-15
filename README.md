
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

Two `.yml` files are included into the repository at `envs/scripts`. For the full functionality (i.e. analysis and plotting), both of them need to be created via

```bash
conda env create -f NanoporeAmpliconGenotyping/envs/nageno.yml
conda env create -f NanoporeAmpliconGenotyping/envs/nageno_plot.yml
conda activate nageno
```

> [!NOTE]
> Alternative to `conda`, `mamba` or `micromamba` can also be used for the creation of the environments which will be much faster. Since `nageno` ultimately uses the environments for some parts of the analysis, make sure that the command `conda activate nageno` works. Alternatively, make sure that you are setting `--manager` to your respective dependency manager (e.g. `micromamba`) while using the pipeline. 

Further, the somatic variant caller, **ClairS-TO**, and its models need to be installed manually, as explained [here](https://github.com/HKU-BAL/ClairS-TO). 

> [!WARNING]
> `clairs-to` searches for the models at `echo ${CONDA_PREFIX}/bin`. This unfortunately can not be changed easily and thus you need to make sure that `clairs-to_models`, `clairs-to_databases`, and `clairs-to_cna_data` exist in the bin-folder of the `nageno` environment. You can prevent this extra step by, as described above, activating the `nageno` environment first and then proceed with the manual `clairs-to` installation.

<details>
<summary>Condensed relevant information about the manual installation of ClairS-TO (click to expand)</summary>

```bash
# in case of a timeout error (Download error (28) Timeout was reached) try modifying timeout settings (works exactly like this only for mamba and conda)
#conda config --set remote_connect_timeout_secs 30
#conda config --set remote_read_timeout_secs 30 

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

Remember to deactivate the nageno env before using NAGeno.

```bash
conda deactivate
```


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
  --ref                REF                  Reference genome file (.fa file needed, .fai files needs to be present too)
  --bed                BED                  BED file for reference
  --txfile             TXFILE               File for visualization

Optional arguments:
  --manager            MANAGER              Package manager used to activate environments (default: conda)
  --threads            THREADS              Number of cores to use (default: 1)
  --min-q              MIN_Q                Minimum base quality (default: 30)
  --max-u              MAX_U                Percentage of bases allowed below MIN_Q (default: 10)
  --mapq               MAPQ                 Minimum mapping quality (default: 50)
  --analysis-dir       DIR                  Directory for output (default: ./analysis)
  --ext                EXT                  Sample name extension (default: SQK-RBK114-24_barcode)
  --clairs-to-path     CLAIR_PATH           Absolute path to 'run_clairs_to' - depends on where ClairS-TO was installed. (default: run_clairs_to)
  --clairs-to-model    CLAIR_MODEL          Clairs-to model (default: ont_r10_dorado_sup_5khz)
  --snpeff-ref         SNPEFF_REF           SNPeff reference genome - should always be the same as the one used for alignment (default: GRCh38.p14)

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
  --min-q              MIN_Q                Minimum base quality (default: 30)
  --max-u              MAX_U                Percentage of bases allowed below MIN_Q (default: 10)
  --mapq               MAPQ                 Minimum mapping quality (default: 50)
  --analysis-dir       DIR                  Directory for output (default: ./analysis)
  --clairs-to-model    CLAIR_MODEL          Clairs-to model (default: ont_r10_dorado_sup_5khz)

```

## Tutorial

Using the exemplary test data in `tutorial`, the correct setup can be confirmed and exemplary output can be generated:

```bash
nageno analysis \
  --dir tutorial/test_data/fastq \
  --anno tutorial/Src/barcode_assignment.tsv \
  --ref /path/to/ref/genome/hg38.fa \
  --bed tutorial/Src/geno_panel_v4.1.bed \
  --txfile tutorial/Src/tx.tsv \
  --analysis-dir tutorial/analysis \
  --threads 20 \
  --clairs-to-path /path/to/run_clairs_to
```

<details>
<summary>Potential installation errors:</summary>

- `[ERROR] file .../envs/nageno/bin/clairs-to_models/ont_r10_dorado_sup_5khz/pileup_affirmative.pkl not found`: Make sure that `clairs-to_models`, `clairs-to_databases`, and `clairs-to_cna_data` exist in the bin-folder of the `nageno` environment. => The best way to ensure that is by installing ClairS-TO while the nageno env is activated.
- `[ERROR] while connecting to https://snpeff.blob.corewindows.net/databases/v5_2snpEff_v5_2[refGenomeVersion].zip`: SnpEff usually downloads the required databases automatically, however, every few years they re-structure which can lead to issues. Try a manual download within the nageno env at `.../mamba/envs/nageno/share/snpeff-5.2-1/` via `snpEff -download [refGenomeVersion]` or use another database. All databases can be viewed there via `snpEff databases`. The annotation database should always match the database previously used for annotation and variant calling. You can read more on that issue [here](https://www.biostars.org/p/296349/). 
- ...

</details>


The `nageno plot` subfunction results in the creation of various different visualisations for the `nageno analysis` output. This is supposed to be used as a quick and comprehensive overview about the genotypes of your samples.

```bash

nageno plot \
  --dir tutorial/test_data/fastq \
  --anno tutorial/Src/barcode_assignment.tsv \
  --ref /path/to/ref/genome/hg38.fa \
  --bed tutorial/Src/geno_panel_v4.1.bed \
  --txfile tutorial/Src/tx.tsv \
  --analysis-dir tutorial/analysis \
  --threads 20 \
  --clairs-to-path /path/to/run_clairs_to
```

> [!TIP]
> `nageno plot` needs less arguments than `nageno analysis`. Since additional arguments are ignored, the quickest way to use the plotting functionality on your results is by replacing the `analysis` with the `plot` subcommand and re-run.

Here is an overview of the generated output files for the provided test data.
Additionally, per sample vcf files, more elaborate vcf collection files for all samples, filtered fastq and filtered bam files will be saved along with log files and htmls files generated by fastplong and SnpEff.

#### Table 1 – SNV genotyping results (all SNVs): `SNV_genotyping_results.tsv`

#### Table 2 – SNV genotyping results (protein-coding SNVs): `prot_coding_SNV_genotyping_results.tsv`

#### Table 3 – Summary of depth statistics: `Summary_depth_stats.tsv`



#### Plot 1 – Allele frequency plot (protein-coding SNVs): `prot_coding_SNV_genotyping_results.svg`

#### Plot 2 – Allele frequency plot (all SNVs): `SNV_genotyping_results.svg`

#### Plots 3-x – Per sample, per gene depth plots: `[GENE]_by_sample_[FILTER].svg` and `[GENE]_[FILTER].svg`


## Citation and Contribution
BioRXive link / doi

## License

This project is licensed under the [Apache License 2.0](LICENSE).
[![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg)](LICENSE)

NOTES: Outdated?
- currently needs to be started NanoporeAmpliconGenotyping dir to work - otherwise it doesn't find the scripts (ERROR: /tmp/tmpyzxjvtjl: line 3: ./scripts/genoSuperscript.sh: No such file or directory)
- how about introducing a starting step like --start-at var-calling and internal options to start at any step in the pipeline?
###############################################################################################
