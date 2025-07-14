# genotype UM based on np data
# plot port coding muts w GQ>=10 and export filtered cols of all data (genotyping_results.tsv) and only prot coding muts w GQ>=10 (genotyping_prot_coding_results.tsv) 

suppressPackageStartupMessages({
    library(readr)
    library(tidyr)
    library(dplyr)
    library(ggplot2)
    library(viridis)
    library(RColorBrewer)
    library(svglite)
    library(hrbrthemes)
})

############

# colors:


col<- hcl.colors(n=7, palette = "viridis")
col
#Farben definieren
col <- hcl.colors(n=20, palette = "viridis")
col

#Farben
colortest <- data.frame("rows"=rep(1,20))
colortest
# barplot(colortest[,1], col = col, names.arg = c(1:20))


# barplot(1:3, col = col[c(2, 10, 17)], names.arg = c(3, 13, 17))
col_3 = col[c(2, 10, 17)]
col_4 = col[c(2, 10, 17, 20)]

###################################


#  load data
args <- commandArgs(trailingOnly = TRUE)

# first arg is dir, second is tsv file
dir = args[1]
#dir = "/home/vera/gueseer/Projects/cancerna/genotyping/np_amplicon_geno/analysis/geno_all_sr/no_trim/sup/"


data_file = args[2]
#data_file = "q90_Q30_MAPQ50_ss_vcf_collection_geno2-3.tsv"
mod = args[3]
#mod = "q90_Q30_MAPQ50_ss"


data = data.frame(read_tsv(file = paste0(dir, data_file), col_names = T))
print(data_file)

print("Loaded data.")

data$SAMPLE_STRIPPED = gsub(pattern ="_.*", "", x = data$SAMPLE)
data$seqed_gene = data$GENE
data[data$seqed_gene != "gene", "seqed_gene"] = "black"
data[data$seqed_gene == "gene", "seqed_gene"] = "red"


data$CHROM = factor(data$CHROM)
data$GENE = factor(data$GENE)
data$SAMPLE = factor(data$SAMPLE)
data$FILTER = factor(data$FILTER)
data$Annotation = factor(data$Annotation)
data$Annotation_Impact = factor(data$Annotation_Impact)
data$Feature_Type = factor(data$Feature_Type)
data$Transcript_BioType = factor(data$Transcript_BioType)
data$SAMPLE_STRIPPED = factor(data$SAMPLE_STRIPPED)
data$seqed_gene = factor(data$seqed_gene)
data$mutGeneID.c = factor(paste0(data$HGVS.c, "_", data$GENE))
data$mutGeneID.p = factor(with(data, ifelse(is.na(HGVS.p), NA, paste0(HGVS.p, "_", GENE))))


head(data)

# re-arrange cols
data = data[,c("SAMPLE_STRIPPED", "GENE", "CHROM", "HGVS.p", "HGVS.c", "AF", "POS", "REF", "ALT", "GQ", "DP", "FILTER", "GT", "AD", "QUAL", "Annotation", "Annotation_Impact", "Feature_ID", "mutGeneID.p", "mutGeneID.c")]

# subset table for plotting
results = subset(data, GQ >= 10)

# results for prot-coding only
res_prot_cod = subset(results, (!is.na(HGVS.p) & GQ>=10))

# fill missing values with AF to make plot better - assumes that whole amplicon is covered for each sample - potentially include a testing step here, but previous QC should cover thsi already
# Complete only missing SAMPLE_STRIPPED values for existing (HGVS.p, GENE)
results_plot <- results %>%
  complete(nesting(HGVS.c, GENE), SAMPLE_STRIPPED, fill = list(AF = 0))

res_prot_cod_plot <- res_prot_cod %>%
  complete(nesting(HGVS.p, GENE), SAMPLE_STRIPPED, fill = list(AF = 0))

#############
### TO DO ###
#############
# consider sample size and number of SNVs (w adn w/o intronic) for plot measure ments and facet_wrap
# make one pc and one CDS plot
# clean code

# test different sample sizes
# pc_test = res_prot_cod_plot
# cds_test = results_plot

# res_prot_cod_plot = subset(pc_test, SAMPLE_STRIPPED == "M15")
# results_plot = subset(cds_test, SAMPLE_STRIPPED == "M15")
# unique(results_plot$mutGeneID.p)

print("#################################################################")
print("Start plotting.")

pc_sample_count = length(unique(res_prot_cod_plot$SAMPLE_STRIPPED))
cds_sample_count = length(unique(results_plot$SAMPLE_STRIPPED))

pc_snv_count = sum(!is.na(unique(res_prot_cod_plot$mutGeneID.p)))
cds_snv_count = sum(!is.na(unique(results_plot$mutGeneID.c)))

pc_combined_counter = pc_sample_count * pc_snv_count
cds_combined_counter = cds_sample_count * cds_snv_count

plot_pc = ggplot(res_prot_cod_plot, aes(x = SAMPLE_STRIPPED, y = AF, fill = SAMPLE_STRIPPED))+
  geom_bar(stat = "identity", position = "dodge")+
  theme_ipsum()+
  ggtitle("Genotyping results - protein coding SNVs")+
  geom_hline(yintercept = 0.5)+
  ylab("Allele Frequency")+
  xlab("Sample")+
  scale_fill_manual(values = viridis(pc_sample_count))+
  scale_y_continuous(breaks = scales::pretty_breaks(n=3))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.position = "none",
        axis.title.x = element_text(hjust = 0.5, size = 14, margin = margin(t = 15)),
        axis.title.y = element_text(hjust = 0.5, size = 14, margin = margin(r = 15)))+
  facet_wrap(~GENE+HGVS.p, scale = "free_x", ncol = round(sqrt(sqrt(pc_combined_counter))+1.7-(pc_sample_count/pc_snv_count)))

plot_cds = ggplot(results_plot, aes(x = SAMPLE_STRIPPED, y = AF, fill = SAMPLE_STRIPPED))+
  geom_bar(stat = "identity", position = "dodge")+
  theme_ipsum()+
  ggtitle("Genotyping results - all SNVs")+
  geom_hline(yintercept = 0.5)+
  ylab("Allele Frequency")+
  xlab("Sample")+
  scale_fill_manual(values = viridis(cds_sample_count))+
  scale_y_continuous(breaks = scales::pretty_breaks(n=3))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.position = "none",
        axis.title.x = element_text(hjust = 0.5, size = 14, margin = margin(t = 15)),
        axis.title.y = element_text(hjust = 0.5, size = 14, margin = margin(r = 15)))+
  facet_wrap(~GENE+HGVS.c, scale = "free_x", ncol = round(sqrt(sqrt(cds_combined_counter))+1.7-(cds_sample_count/cds_snv_count)))

print("Save plots.")

# save result as tsv
write_tsv(data, file = paste0(dir, mod, "-SNV_genotyping_results.tsv"))
write_tsv(res_prot_cod, file = paste0(dir, mod, "-prot_coding_SNV_genotyping_results.tsv"))

# save plots
svglite(filename = paste0(dir, mod, "-prot_coding_SNV_genotyping_results.svg"),
        width = (round(sqrt(sqrt(pc_combined_counter))+2.5-(pc_sample_count/pc_snv_count)) * (pmax(pc_sample_count, 3) * 0.6) / sqrt(sqrt(pc_combined_counter)))+4,
        height = pmax(ceiling(pc_snv_count / round(sqrt(sqrt(pc_combined_counter))+1.7-(pc_sample_count/pc_snv_count)))*2.5, 4))
print(plot_pc) # print ensures the plot is actually printed, otherwise timeout before saving
dev.off()

svglite(filename = paste0(dir, mod, "-SNV_genotyping_results.svg"),
        width = (round(sqrt(sqrt(cds_combined_counter))+2.5-(cds_sample_count/cds_snv_count)) * (pmax(cds_sample_count, 3) * 0.6) / sqrt(sqrt(cds_combined_counter)))+4,
        height = ceiling(cds_snv_count / round(sqrt(sqrt(cds_combined_counter))+1.7-(cds_sample_count/cds_snv_count)))*2.5)
print(plot_cds) # print ensures the plot is actually printed, otherwise timeout before saving
dev.off()

print(paste0("Output saved in", dir))
