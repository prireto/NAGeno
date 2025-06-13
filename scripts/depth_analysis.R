# visualize seqdepth

# install.packages("fuzzyjoin", repos = "https://cloud.r-project.org/")

library(ggplot2)
library(hrbrthemes)
library(readr) # import tsv files
library(Formula)
library(ggpubr)
library(viridis) # pretty colors
library(ggrepel) # better way of preventing label and text overlaps in ggplot2
library(scales)
library(svglite) # save plots as svgs
library(rlang) # for interpreting strings as symbols ( enables pasted strings to be var placeholders)
library(dplyr)
library(fuzzyjoin)

args <- commandArgs(trailingOnly = TRUE)

files = args[1]
print(files)
# files="/home/vera/gueseer/Pipelines/NanoporeAmpliconGenotyping/test_data/filtered_bam_sr/depth/"

files_specifier = args[2]
print(files_specifier)
# files_specifier = "SQK-RBK114-24_barcode"

files_mod = args[3]
print(files_mod)
# files_mod = "_q90_Q30_MAPQ50"
# files_mod = ""

plot_dir = paste0(files, "plots/")

# load bc annotation
anno_dir = args[4]
print(anno_dir)
# anno_dir = "/home/vera/gueseer/Pipelines/NanoporeAmpliconGenotyping/test_data/anno_trial.tsv"
bc_anno = data.frame(read_tsv(file = anno_dir, col_names = c("sample", "BC")))
# Convert the 'samples' column in anno to a factor with correctly ordered levels
bc_anno$sample = factor(bc_anno$sample, levels = bc_anno$sample[order(as.numeric(sub("S", "", bc_anno$sample)))])
head(bc_anno)

# samples = c("Mel202", "MH081", "MH080", "MH079", "MH037", "MH025", "MH023", "MH022", "MH006", "MH002")

mut_list_dir = args[5]
print(mut_list_dir)
# mut_list_dir = "/home/vera/gueseer/Pipelines/NanoporeAmpliconGenotyping/test_data/mut_list.tsv"
muts = data.frame(read_tsv(file = paste0(mut_list_dir), col_names = T))
muts


bed_dir = args[6]
print(bed_dir)
# bed_dir = "/home/vera/gueseer/Src/geno_panel_v4.1.bed"
bed = data.frame(read_tsv(file = bed_dir, col_names = c("chr", "start", "stop", "gene")))
bed

# read data, add sample col and combine all into one df
data <- list()
for (samp in bc_anno[["sample"]]) {
  filename = paste0(files_specifier, bc_anno[bc_anno$sample == samp, "BC"], "_", samp, files_mod, ".sorted.bam.depth.tsv")
  data <- rbind(data, cbind(data.frame(read_tsv(file = paste0(files, filename), col_names = T)), sample = samp))
}
head(data)
unique(data$sample)


# Perform a range-based join to assign gene names based on bed - this might take a wile
data <- fuzzy_inner_join(
  data,
  bed,
  by = c(
    "contig" = "chr",               # exact match on contig/chr
    "position" = "start",           # lower bound of range
    "position" = "stop"             # upper bound of range
  ),
  match_fun = list(`==`, `>=`, `<=`)  # contig must match, position must be within [start, stop]
) %>%
  select(-chr, -start, -stop)  # clean up
head(data)

genes = unique(data$gene)

# ggplot(dataQ, aes(x=position, y= depth, fill = sample))+
#   geom_line()+
#   theme_ipsum()+
#   facet_wrap(~gene, scales = "free_x")


#############
### TO DO ###
#############

# omit all by sample plots? or create them optionally?
# scale export size based on sample size for by sample plots
# set y axis scale based on highest occurring depth
# adjust plot size for cumulative plot based on legend - one common plot?

# test plotting with many samples:

data_test = data
for (i in 1:25) {
  test = data
  test[test$sample == "S3", "sample"] = paste0("T", i)
  test[test$sample == "S8", "sample"] = paste0("T", i+50)
  data_test = rbind(data_test, test)
  # print(head(test))
}
unique(data_test$sample)
length(unique(data_test$sample))


#############
### PLOTS ###
#############

# check depth for each gene by sample
for (gene in genes) {
  plot_by_sample = ggplot(subset(data_test, gene == gene), aes(x = position, y = depth, color = sample))+
    geom_line()+
    theme_ipsum()+
    ggtitle(paste0(gene, " amplicon - per base depth"))+
    scale_color_manual(values = viridis(length(unique(data_test$sample))))+
    scale_y_continuous(trans='log10', limits = c(1, 100000), labels = comma)+
    facet_wrap(~sample, ncol = round(sqrt(length(unique(data_test$sample)))))+
    geom_hline(yintercept = 10, color = "red", linetype = "dashed")+
    xlim(c(bed[bed$gene == gene, "start"], bed[bed$gene == gene, "stop"]))+
    theme(axis.text.x = element_text(angle = 90), legend.position = "none",
          axis.title.x = element_text(hjust = 0.5, size = 14),
          axis.title.y = element_text(hjust = 0.5, size = 14)
          )
  print(plot_by_sample)
  
  # save plot
  file = paste0(plot_dir, "amplicon_geno_", gene, "_by_sample", files_mod)
  svglite(filename = paste0(file, ".svg"), width = round(sqrt(length(unique(data_test$sample))))+4, height = round(sqrt(length(unique(data_test$sample))))+4)
  print(plot_by_sample)
  dev.off()
}

# check depth for each gene in one plot
plot_collection = ggplot()
for (gene in genes) {
  # generate plot
  plot_by_gene = ggplot(subset(data_test, gene == gene), aes(x = position, y = depth, color = sample))+
    geom_line()+
    theme_ipsum()+
    ggtitle(paste0(gene, " amplicon - per base depth"))+
    scale_color_manual(values = viridis(length(unique(data_test$sample))),
                       guide = guide_legend(ncol = pmax(round(sqrt(length(unique(data_test$sample)))) - round(sqrt(length(unique(data_test$sample)))/2), 5))
                       )+
    scale_y_continuous(trans='log10', limits = c(1, 100000), labels = comma)+
    geom_hline(yintercept = 10, color = "red", linetype = "dashed")+
    xlim(c(bed[bed$gene == gene, "start"], bed[bed$gene == gene, "stop"]))+
    theme(axis.text.x = element_text(angle = 90), legend.position = "bottom",
          axis.title.x = element_text(hjust = 0.5, size = 14),
          axis.title.y = element_text(hjust = 0.5, size = 14)
          )
  # print(plot_by_gene)
  
  # save plot
  file = paste0(plot_dir, "amplicon_geno_", gene, files_mod)
  svglite(filename = paste0(file, ".svg"), width = 5, height = 5 +(round(sqrt(length(unique(data_test$sample))))/3*2))
  print(plot_by_gene)
  dev.off()
}

###################
### DEPTH STATS ###
###################
stats = data.frame()

for (gene in genes) {
  chr = bed[bed$gene == gene, "chr"]
  min = min(subset(data, contig == chr & position > bed[bed$chr == chr, "start"] & depth > 1)$position)
  max = max(subset(data, contig == chr & position < bed[bed$chr == chr, "stop"] & depth > 1)$position)
  
  stats = rbind(stats, subset(data, contig == chr & position > min & position < max))
}


stats_summary <- stats %>%
  select(sample, contig, depth, gene) %>%
  group_by(sample, contig, gene) %>%
  summarise(
    median = median(depth, na.rm = TRUE),
    mean = round(mean(depth, na.rm = TRUE), 2),
    # min = min(depth, na.rm = TRUE),
    # max = max(depth, na.rm = TRUE),
    .groups = "drop"
  )
head(stats_summary)


write_tsv(stats_summary, file = paste0(plot_dir, files_mod, "-depth_stats.tsv"))
