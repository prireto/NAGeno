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
    library(scales)
    library(showtext)
})

############

# add narrow font for plotting
font_add_google("Roboto Condensed", "roboto_condensed")
showtext_auto()

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
# dir = "/home/vera/gueseer/Pipelines/NanoporeAmpliconGenotyping/tutorial/analysis/ClairS-TO/"

data_file = args[2]
# data_file = "q90_Q30_MAPQ50_ont_r10_dorado_sup_5khz_ss_snv_vcf_collection.tsv"

out_dir = args[3]
# out_dir = "/home/vera/gueseer/Pipelines/NanoporeAmpliconGenotyping/tutorial/analysis/output/"


data = data.frame(read_tsv(file = paste0(dir, data_file), col_names = T))
print(data_file)

print("Loaded data.")

data$SAMPLE = gsub(pattern ="_.*", "", x = data$SAMPLE_LONG)
data$seqed_gene = data$GENE
data[data$seqed_gene != "gene", "seqed_gene"] = "black"
data[data$seqed_gene == "gene", "seqed_gene"] = "red"


data$CHROM = factor(data$CHROM)
data$GENE = factor(data$GENE)
data$SAMPLE_LONG = factor(data$SAMPLE_LONG)
data$FILTER = factor(data$FILTER)
data$Annotation = factor(data$Annotation)
data$Annotation_Impact = factor(data$Annotation_Impact)
data$Feature_Type = factor(data$Feature_Type)
data$Transcript_BioType = factor(data$Transcript_BioType)
data$SAMPLE = factor(data$SAMPLE)
data$seqed_gene = factor(data$seqed_gene)
data$mutGeneID.c = factor(paste0(data$HGVS.c, "_", data$GENE))
data$mutGeneID.p = factor(with(data, ifelse(is.na(HGVS.p), NA, paste0(HGVS.p, "_", GENE))))


head(data)

# re-arrange cols
data = data[,c("SAMPLE", "GENE", "CHROM", "HGVS.p", "HGVS.c", "AF", "POS", "REF", "ALT", "GQ", "DP", "FILTER", "GT", "AD", "QUAL", "Annotation", "Annotation_Impact", "Feature_ID", "mutGeneID.p", "mutGeneID.c")]

# subset table for plotting
results = subset(data, GQ >= 10)

# results for prot-coding only
res_prot_cod = subset(results, (!is.na(HGVS.p) & GQ>=10))

# fill missing values with AF to make plot better - assumes that whole amplicon is covered for each sample - potentially include a testing step here, but previous QC should cover this already
# Complete only missing SAMPLE values for existing (HGVS.p, GENE)
results_plot <- results %>%
  complete(nesting(HGVS.c, GENE), SAMPLE, fill = list(AF = 0))

res_prot_cod_plot <- res_prot_cod %>%
  complete(nesting(HGVS.p, GENE), SAMPLE, fill = list(AF = 0))


###########################################################
# test cases

# data_save = res_prot_cod_plot
# data_s3 = subset(data_save, SAMPLE == "S3")
# data_G11 = subset(data_save, GENE == "GNA11" & HGVS.p == "p.Thr257Thr")
# ####
# 
# 
# res_prot_cod_plot = data_save
# for (x in 1:48) {
#   data_tmp = transform(data_s3, SAMPLE = paste0("E", x))
#   res_prot_cod_plot = rbind(res_prot_cod_plot, data_tmp)
# }
# unique(res_prot_cod_plot$SAMPLE)
# 
# for (x in 1:24) {
#   data_tmp = transform(data_G11, GENE = paste0("GNA", x))
#   res_prot_cod_plot = rbind(res_prot_cod_plot, data_tmp)
# }
# unique(res_prot_cod_plot$SAMPLE)
# unique(res_prot_cod_plot$GENE)
####################################################

print("#################################################################")
print("Start plotting.")

pc_sample_count = length(unique(res_prot_cod_plot$SAMPLE))
cds_sample_count = length(unique(results_plot$SAMPLE))

pc_snv_count = sum(!is.na(unique(res_prot_cod_plot$mutGeneID.p)))
cds_snv_count = sum(!is.na(unique(results_plot$mutGeneID.c)))

pc_combined_counter = pc_sample_count * pc_snv_count
cds_combined_counter = cds_sample_count * cds_snv_count
###########################################################################
# calc ncol nrow and plot sizes for protein-coding plot
n_facets <- length(unique(paste(res_prot_cod_plot$GENE, res_prot_cod_plot$HGVS.p)))
samples_per_facet <- pc_sample_count


# Layout params
bar_width_per_sample <- 0.28
min_facet_width <- 2.8
facet_width <- max(samples_per_facet * bar_width_per_sample, min_facet_width)

# Plot layout
target_total_width <- 16     # inches
max_cols <- 6                # prevent too many columns
min_cols <- 2

ncol <- max(min_cols, min(max_cols, floor(target_total_width / facet_width)))
nrow <- ceiling(n_facets / ncol)

# Final dimensions
facet_height <- 2.2
title_buffer <- 1.2
plot_width <- ncol * facet_width
plot_height <- nrow * facet_height + title_buffer

###########################################################################

plot_pc = ggplot(res_prot_cod_plot, aes(x = SAMPLE, y = AF, fill = SAMPLE))+
  geom_bar(stat = "identity", position = "dodge")+
  theme_ipsum(base_family = "roboto_condensed")+
  ggtitle("Genotyping results - protein-coding SNVs")+
  geom_hline(yintercept = 0.5)+
  ylab("Allele Frequency")+
  xlab("Sample")+
  scale_fill_manual(values = viridis(pc_sample_count))+
  scale_y_continuous(breaks = scales::pretty_breaks(n=3))+
  facet_wrap(~GENE+HGVS.p, scale = "free_x", ncol = ncol)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.position = "none",
        axis.title.x = element_text(hjust = 0.5, size = 14, margin = margin(t = 15)),
        axis.title.y = element_text(hjust = 0.5, size = 14, margin = margin(r = 15)))
# ncol = round(sqrt(sqrt(pc_combined_counter))+1.7-(pc_sample_count/pc_snv_count))

# save plot
svglite(filename = paste0(out_dir, "Prot_coding_SNV_genotyping_results.svg"),
        width = plot_width,
        height = plot_height)
print(plot_pc) # print ensures the plot is actually printed, otherwise timeout before saving
dev.off()



###########################################################################
# calc ncol nrow and plot sizes for cds plot
n_facets <- length(unique(paste(results_plot$GENE, results_plot$HGVS.c)))
samples_per_facet <- pc_sample_count


# Layout params
bar_width_per_sample <- 0.28
min_facet_width <- 2.8
facet_width <- max(samples_per_facet * bar_width_per_sample, min_facet_width)

# Plot layout
target_total_width <- 16     # inches
max_cols <- 6                # prevent too many columns
min_cols <- 2

ncol <- max(min_cols, min(max_cols, floor(target_total_width / facet_width)))
nrow <- ceiling(n_facets / ncol)

# Final dimensions
facet_height <- 2.2
title_buffer <- 1.2
plot_width <- ncol * facet_width
plot_height <- nrow * facet_height + title_buffer

###########################################################################

plot_cds = ggplot(results_plot, aes(x = SAMPLE, y = AF, fill = SAMPLE))+
  geom_bar(stat = "identity", position = "dodge")+
  theme_ipsum(base_family = "roboto_condensed")+
  ggtitle("Genotyping results - all SNVs")+
  geom_hline(yintercept = 0.5)+
  ylab("Allele Frequency")+
  xlab("Sample")+
  scale_fill_manual(values = viridis(cds_sample_count))+
  scale_y_continuous(breaks = scales::pretty_breaks(n=3))+
  facet_wrap(~GENE+HGVS.c, scale = "free_x", ncol = ncol)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.position = "none",
        axis.title.x = element_text(hjust = 0.5, size = 14, margin = margin(t = 15)),
        axis.title.y = element_text(hjust = 0.5, size = 14, margin = margin(r = 15)))
  
  
# save plot
svglite(filename = paste0(out_dir, "SNV_genotyping_results.svg"),
        width = plot_width,
        height = plot_height)
print(plot_cds) # print ensures the plot is actually printed, otherwise timeout before saving
dev.off()
###########################################################################

print("Save results.")

# save result as tsv
write_tsv(data, file = paste0(out_dir, "SNV_genotyping_results.tsv"))
write_tsv(res_prot_cod, file = paste0(out_dir, "Prot_coding_SNV_genotyping_results.tsv"))


print(paste0("Output saved in ", out_dir))
