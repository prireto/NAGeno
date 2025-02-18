# genotype UM based on np data
# plot port coding muts w GQ>=10 and export filtered cols of all data (genotyping_results.tsv) and only prot coding muts w GQ>=10 (genotyping_prot_coding_results.tsv) 

library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(svglite)
library(hrbrthemes)

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
# dir = "/home/vera/gueseer/Projects/cancerna/genotyping/np_amplicon_geno/analysis/geno1_sr/ClairS-TO/"


data_file = args[2]
# data_file = "2025-01-02_q90_Q30_MAPQ50_ss_vcf_collection_geno2-3.tsv"
mod = args[3]
# mod = "q90_Q30_MAPQ50_ss"


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



# flag if either REF=C and ALT=T or REF=G and ALT=A
data$deamination = FALSE
data[((data$REF == "C" & data$ALT == "T") | (data$REF == "G" & data$ALT =="A")), "deamination"] = TRUE

head(data)

# re-arrange cols
data = data[,c("SAMPLE_STRIPPED", "GENE", "CHROM", "HGVS.p", "HGVS.c", "AF", "POS", "REF", "ALT", "GQ", "DP", "FILTER", "GT", "AD", "QUAL", "Annotation", "Annotation_Impact", "Feature_ID", "deamination")]

# subset table for results
results = subset(data)

# results for prot-coding only
res_prot_cod = subset(results, (!is.na(HGVS.p) & GQ>=10))


#  fix this - SRP
# # complete data for plotting
# # Define all combinations of samples and prot level muts
# all_combinations = expand.grid(
#   SAMPLE_STRIPPED = unique(res_prot_cod$SAMPLE_STRIPPED),
#   HGVS.p = unique(res_prot_cod$HGVS.p)
# )
# 
# # Fill in missing rows with NA for value
# data_ext = res_prot_cod %>%
#   full_join(all_combinations, by = c("SAMPLE_STRIPPED", "HGVS.p", "GENE")) %>%
#   mutate(AF = replace_na(AF, 0)) # Replace NA with 0



# fill missing values with AF to make plot better - assumes that whole amplicon is covered for each sample - potentially include a testing step here, but previous QC should cover thsi already
# Complete only missing SAMPLE_STRIPPED values for existing (HGVS.p, GENE)
res_prot_cod_plot <- res_prot_cod %>%
complete(nesting(HGVS.p, GENE), SAMPLE_STRIPPED, fill = list(AF = 0))


print("#################################################################")
print("Start plotting.")


# Now plot by gene 
plot = ggplot(res_prot_cod, aes(x = HGVS.p, y = AF, fill = SAMPLE_STRIPPED))+
  # geom_bar(position = position_dodge2(preserve = "single"), stat = "identity")+
  geom_bar(stat = "identity", position = "dodge")+
  theme_ipsum()+
  ggtitle("Geno results")+
  geom_hline(yintercept = 0.5)+
  ylab("Allele Frequency")+
  xlab("protein-coding SNV")+
  scale_fill_manual(values = viridis(length(unique(res_prot_cod[,"SAMPLE_STRIPPED"]))))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  facet_wrap(~GENE, axes = "all_x", scales = "free_x")
# plot

plot2 = ggplot(res_prot_cod_plot, aes(x = HGVS.p, y = AF, fill = SAMPLE_STRIPPED))+
  # geom_bar(position = position_dodge2(preserve = "single"), stat = "identity")+
  geom_bar(stat = "identity", position = "dodge")+
  theme_ipsum()+
  ggtitle("Geno results")+
  geom_hline(yintercept = 0.5)+
  ylab("Allele Frequency")+
  xlab("protein-coding SNV")+
  scale_fill_manual(values = viridis(length(unique(res_prot_cod[,"SAMPLE_STRIPPED"]))))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  facet_wrap(~GENE, axes = "all_x", scales = "free_x")
# plot2


plot3 = ggplot(res_prot_cod_plot, aes(x = SAMPLE_STRIPPED, y = AF, fill = SAMPLE_STRIPPED))+
  # geom_bar(position = position_dodge2(preserve = "single"), stat = "identity")+
  geom_bar(stat = "identity", position = "dodge")+
  theme_ipsum()+
  # ggtitle("Geno1 results")+
  geom_hline(yintercept = 0.5)+
  ylab("Allele Frequency")+
  xlab("Sample")+
  scale_fill_manual(values = viridis(length(unique(res_prot_cod[,"SAMPLE_STRIPPED"]))))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.position = "none",
        axis.title.x = element_text(hjust = 0.5, size = 14, margin = margin(t = 15)),
        axis.title.y = element_text(hjust = 0.5, size = 14, margin = margin(r = 15)))+
  facet_wrap(~GENE+HGVS.p, scale = "free_x")
# plot3

print("Save plots.")


# save result as tsv
write_tsv(results, file = paste0(dir, mod, "-genotyping_results.tsv"))
write_tsv(res_prot_cod, file = paste0(dir, mod, "-genotyping_prot_coding_results.tsv"))

# save plots
svglite(filename = paste0(dir, mod, "-genotyping_results.svg"), width = 15, height = 10.5)
print(plot) # print ensures the plot is actually printed, otherwise timeout before saving
dev.off()

svglite(filename = paste0(dir, mod, "-skinny-genotyping_results.svg"), width = 15, height = 10.5)
print(plot2) # print ensures the plot is actually printed, otherwise timeout before saving
dev.off()

svglite(filename = paste0(dir, mod, "-by_mut-genotyping_results.svg"), width = 15, height = 12)
print(plot3) # print ensures the plot is actually printed, otherwise timeout before saving
dev.off()

print(paste0("Output saved in", dir))

# 
# 
# # Original test_data with missing rows for some sample-category combinations
# test_data <- data.frame(
#   sample = c("Sample1", "Sample1", "Sample2", "Sample3"),
#   category = c("A", "B", "A", "C"),
#   value = c(10, 15, 5, 20)
# )
# 
# # Define all combinations of samples and categories
# all_combinations <- expand.grid(
#   sample = unique(test_data$sample),
#   category = unique(test_data$category)
# )
# 
# # Fill in missing rows with NA for value
# test_data_complete <- test_data %>%
#   full_join(all_combinations, by = c("sample", "category")) %>%
#   mutate(value = replace_na(value, 0)) # Replace NA with 0 (or leave as NA if preferred)
# 
# # Grouped barplot
# ggplot(test_data_complete, aes(x = category, y = value, fill = sample)) +
#   geom_bar(stat = "identity", position = position_dodge()) +
#   theme_minimal() +
#   labs(title = "Grouped Barplot with Gaps for Missing Values",
#        x = "Category",
#        y = "Value") +
#   scale_fill_brewer(palette = "Set3")
# 



