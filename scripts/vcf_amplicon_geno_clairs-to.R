#analyse vcf files from np_amplicon_geno
#extract only result lines and relevant tx, also separate genotype/allele freq col into single cols
suppressPackageStartupMessages({
    library(readr)
    library(tidyr)
    library(dplyr)
    library(stringr)
})

# get SAMPLEs
args <- commandArgs(trailingOnly = TRUE)

# input args should be (1) dir w dirs for snv analysis per sample, (2) modifier (e.g. _q90Q20_ss) and then the names of all samples to be processed

dir = args[1]
# dir = "/home/vera/gueseer/Projects/cancerna/genotyping/np_amplicon_geno/analysis/geno3_sr/no_trim/ClairS-TO/"
print(paste0("Working in ", dir))

mod = args[2]
# mod = "q90_Q30_MAPQ50_ss"

txfile = args[3]
# txfile = "/home/vera/gueseer/Src/geno/tx.tsv"

outdir = args[4]
# outdir = "/home/vera/gueseer/Projects/cancerna/genotyping/np_amplicon_geno/analysis/geno3_sr/no_trim/ClairS-TO/MH002_q90_Q30_MAPQ50_ss/"

type = args[5]
# type = "indel"

samps = args[6:length(args)]
# samps = c("MH022", "MH037", "MH079", "MH080", "MH081", "x92.1", "OMM1.5", "Mel285", "Mel270", "M15", "MH024", "MH026", "MH032", "MH091", "MH010", "MH082", "MH028", "MH001", "MH009", "MH031", "MH027", "MH029", "MH083", "MH084")
# samps = "MH002"

# load tx file
tx = data.frame(read_tsv(file = txfile, col_names = T))

# args = c("Mel202", "Mel202_q90_Q20")
# args = ("Mel202Filtered")
print(tx)
print(samps)
print(paste0("start: ", samps[1]))


#setwd("/home/vera/gueseer/Projects/cancerna/genotyping/np_amplicon_geno/analysis/geno1_sr-geno-ref/wf-human-variation")
#set dir
# dir = "/home/vera/gueseer/Projects/cancerna/genotyping/np_amplicon_geno/analysis/geno1_sr-geno-ref/wf-human-variation/"
# dir = "/home/vera/gueseer/Projects/cancerna/genotyping/np_amplicon_geno/analysis/geno1_sr/ClairS-TO/"

# initiate results var
res = data.frame()


summary(samps)
# different processing for indel and snp vcfs bc indels are not annotated
# for snv, ".anno" needs ot be included in the file name
if(type == "snv"){
  for (samp in samps) {
    SAMPLE=paste0(samp, "_", mod)
    # SAMPLE="Mel202_q90_Q20_ss"
    # load data and sep cols by TAB, make sure ALT and REF ar chars, if they only consist a T they would be declared booleans otherwise and have issues with row binding later on
    print(paste0("hello: ", dir, SAMPLE, "/", type, ".anno.vcf_temp"))
    data = data.frame(read_tsv(paste0(dir, SAMPLE, "/", type, ".anno.vcf_temp"), col_names=T,col_types = cols(ALT = col_character(), REF = col_character())))
    data
    
    # rename first col
    colnames(data)[1] = "CHROM"
    
    # mv both last cols before INFO
    data = data[,c(1:7,9,10,8)]
    # sep "SAMPLE" col by ':', take col names by separating FORMAT col by ':'
    description = data[1,"FORMAT"]
    description = strsplit(description, split=":", fixed=T)
    # is a list now
    description = description[[1]]
    data = separate(data, col=SAMPLE, into = description, sep = ":", remove = T, convert = T)
    data$AD = gsub(",", ";", data$AD)
    # delete FORMAT col
    data = data[,!colnames(data) %in% c("FORMAT")]
    
    # take everything up until "ANN" and save in separate col ADD_INFO (no equal number of cols for the previous info, so skip for now)
    data$ADD_INFO = gsub(pattern = ";ANN=.*", replacement = "", x = data$INFO, perl = TRUE)
    data$INFO = gsub(pattern = ".*;ANN=", replacement = "ANN=", x = data$INFO, perl = TRUE)
    
    # sep cols based on ';' into INFO, ANNO and REST, discard REST - warning
    data = separate(data, col="INFO", into = c("ANNO", "REST"), sep = ";", remove = T, convert = F)
    
    data = data[,!colnames(data) %in% c("REST")]
    # rm "ANN="
    data$ANNO = gsub(pattern = "ANN=", replacement = "", x = data$ANNO)
    # sep cols based on ',' - need to set a size - take sth ridiculously large to be safe, will be deleted then anyways
    # first determine how many new cols are needed
    max_fields = max(str_count(data$ANNO, ",") + 1, na.rm = TRUE)
    # generate dynamic col names
    col_names = paste0("ANNO_", seq_len(max_fields))
    data = separate(data, col = ANNO, into = col_names, sep = ",", remove = FALSE)
    # add first col for SAMPLE
    data = cbind(rep(SAMPLE, length(data[,1])), data)
    colnames(data)[1] = "SAMPLE"
    
    # add col for gene
    data = cbind(rep("gene", length(data[,1])), data)
    colnames(data)[1] = "GENE"
    # if one of the genes in tx$gene is also in one of the cols data[,16 or 17:202] => change "gene" to that gene
    # set "|" as word boundaries for genes with overlapping names
    # !! leads to problems in case of genomic locus overlaps !!
    
    # PROBLEM: in some vcf files there is an additional PL col in some there isn't (maybe based on whether I pre-filter for specific Q scores?)
    # solution: set col_ANNO to make sure it always takes the right cols
    # NEW PROBLEM: need to insert blank PL col in case there is none to match col number - just don't mix pre- and non-pre-filtered files...
    
    col_ANNO = which(colnames(data) %in% "ANNO")
    
    for (i in tx$gene) {
      data[grepl(paste0("\\|", i, "\\|"), data[,col_ANNO]),"GENE"] = i
    }
    
    # rm all post ANNO cols except the one that matches the tx#
    # in all rows of the same gene, search the col containing the right tx and replace ANNO w it - rm all other cols after 16
    # needs for loop to run through data bc not all rows for one gene necessarily have the same number of cols
    for (i in 1:length(tx[,1])) {
      for (j in 1:length(data[,1])) {
        if(data[j, "GENE"]==tx[i, "gene"]){
          data[j, "ANNO"] = data[j, grep(tx[i,"tx"], data[j,-(1:col_ANNO)])+col_ANNO]
        }
      }
    }
    data = data[,1:col_ANNO]
    
    # sep ANNO col by "|" - need to excape | => "\\|"
    anno = c("Allele", "Annotation", "Annotation_Impact", "Gene_Name", "Gene_ID", "Feature_Type", "Feature_ID",
             "Transcript_BioType", "Rank", "HGVS.c", "HGVS.p", "cDNA.pos/cDNA.length", "CDS.pos/CDS.length", "AA.pos/AA.length", "Distance", "ERRORS/WARNINGS/INFO")
    data = separate(data, col=ANNO, into = anno, sep = "\\|", remove = T, convert = F)
    
    # add data of one sample to overall data
    res = dplyr::bind_rows(res, data)
    res$GQ_PASS = ifelse(res$GQ>=10, "PASS", "FAIL")
    print(res)
  }
} else if(type == "indel"){
  for (samp in samps) {
    SAMPLE=paste0(samp, "_", mod)
    # SAMPLE="Mel202_q90_Q20_ss"
    # load data and sep cols by TAB, make sure ALT and REF are chars, if they only consist a T they would be declared booleans otherwise and have issues with row binding later on
    # should also work with Strings (for indels)
    print(paste0("hello: ", dir, SAMPLE, "/", type, ".vcf_temp"))
    data = data.frame(read_tsv(paste0(dir, SAMPLE, "/", type, ".vcf_temp"), col_names=T, col_types = cols(ALT = col_character(), REF = col_character())))
    data
    
    # rename first col
    colnames(data)[1] = "CHROM"
    
    # delete INFO col
    data = data[,!colnames(data) %in% c("INFO")]
    
    # sep "SAMPLE" col by ':', take col names by separating FORMAT col by ':'
    description = data[1,"FORMAT"]
    description = strsplit(description, split=":", fixed=T)
    # is a list now
    description = description[[1]]
    data = separate(data, col=SAMPLE, into = description, sep = ":", remove = T, convert = T)
    data$AD = gsub(",", ";", data$AD)
    
    # delete base count cols
    data = data[,!colnames(data) %in% c("AU", "CU", "GU", "TU")]
    
    # delete FORMAT col
    data = data[,!colnames(data) %in% c("FORMAT")]
    # add first col for SAMPLE
    data = cbind(rep(SAMPLE, length(data[,1])), data)
    colnames(data)[1] = "SAMPLE"
    # add data of one sample to overall data, if non-empty
    if (nrow(data) > 0) {
      res = dplyr::bind_rows(res, data)
    }

    res$GQ_PASS = ifelse(res$GQ>=10, "PASS", "FAIL")
    print(res)
  }
}




# save result as tsv
print(paste0("Output saved in", outdir))
write_tsv(res, file = paste0(outdir, mod, "_", type, "_vcf_collection.tsv"))

