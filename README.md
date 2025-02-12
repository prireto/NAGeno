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





TODO for pipeline
- in general is it possible to make things as options w default settings?
    - genome ref (for alignment, vcf annotation etc.)
    - kmer and window size for minimap
    - max number of threads to use
    - bed file
    - basically important options for clairs-to
- path to tx.tsv needs to be parsed for vcf extraction - no reasonable default here




Notes
- Pipeline needs to run in clirs-to env
- SnpEff needs newer Java version than we have installed
- SnpEff is installed in /home/vera/gueseer/App/tools/snpEff/
