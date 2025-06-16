# visualize seqdepth

suppressPackageStartupMessages({
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
})

args <- commandArgs(trailingOnly = TRUE)
print("Input used for depth plotting:")

files = args[1]
print(files)

files_specifier = args[2]
print(files_specifier)

files_mod = args[3]
print(files_mod)

# load bc annotation
anno_dir = args[4]
print(anno_dir)
bc_anno = data.frame(read_tsv(file = anno_dir, col_names = c("sample", "BC")))

# Convert the 'samples' column in anno to a factor with correctly ordered levels
bc_anno$sample = factor(bc_anno$sample, levels = bc_anno$sample[order(as.numeric(sub("S", "", bc_anno$sample)))])

mut_list_dir = args[5] # currently not used
print(mut_list_dir)

if (mut_list_dir == "NA") {
  mut_list_dir = NULL
} else {
  muts = data.frame(read_tsv(file = paste0(mut_list_dir), col_names = T))
  muts
}

bed_dir = args[6]
print(bed_dir)
bed = data.frame(read_tsv(file = bed_dir, col_names = c("chr", "start", "stop", "gene")))

plot_dir = args[7]
print(plot_dir)

# if any of the input files are empty, stop execution
if (nrow(bc_anno) == 0 || nrow(bed) == 0 || length(files) == 0) {
  stop("One or more input files are empty. Please check the file paths.")
}

# read data, add sample col and combine all into one df
data <- list()
for (samp in bc_anno[["sample"]]) {
  filename = paste0(files_specifier, bc_anno[bc_anno$sample == samp, "BC"], "_", samp, files_mod, ".sorted.bam.depth.tsv")
  data <- rbind(data, cbind(data.frame(read_tsv(file = paste0(files, filename), col_names = T)), sample = samp))
}
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

genes = unique(data$gene)

#############
### PLOTS ###
#############
print("Start depth plotting...")

head(data)
head(genes)

# check depth for each gene by sample
for (gene in genes) {
  plot_by_sample = ggplot(subset(data, gene == gene), aes(x = position, y = depth, color = sample))+
    geom_line()+
    theme_ipsum()+
    ggtitle(paste0(gene, " amplicon - per base depth"))+
    scale_color_manual(values = viridis(length(unique(data$sample))))+
    scale_y_continuous(trans='log10', limits = c(1, 100000), labels = comma)+
    facet_wrap(~sample, ncol = round(sqrt(length(unique(data$sample)))))+
    geom_hline(yintercept = 10, color = "red", linetype = "dashed")+
    xlim(c(bed[bed$gene == gene, "start"], bed[bed$gene == gene, "stop"]))+
    theme(axis.text.x = element_text(angle = 90), legend.position = "none",
          axis.title.x = element_text(hjust = 0.5, size = 14),
          axis.title.y = element_text(hjust = 0.5, size = 14)
          )
  
  # save plot
  file = paste0(plot_dir, "amplicon_geno_", gene, "_by_sample", files_mod)
  svglite(filename = paste0(file, ".svg"), width = round(sqrt(length(unique(data$sample))))+4, height = round(sqrt(length(unique(data$sample))))+4)
  print(plot_by_sample)
  dev.off()
}

# check depth for each gene in one plot
plot_collection = ggplot()

for (gene in genes) {
  # generate plot
  plot_by_gene = ggplot(subset(data, gene == gene), aes(x = position, y = depth, color = sample))+
    geom_line()+
    theme_ipsum()+
    ggtitle(paste0(gene, " amplicon - per base depth"))+
    scale_color_manual(values = viridis(length(unique(data$sample))),
                       guide = guide_legend(ncol = pmax(round(sqrt(length(unique(data$sample)))) - round(sqrt(length(unique(data$sample)))/2), 5))
                       )+
    scale_y_continuous(trans='log10', limits = c(1, 100000), labels = comma)+
    geom_hline(yintercept = 10, color = "red", linetype = "dashed")+
    xlim(c(bed[bed$gene == gene, "start"], bed[bed$gene == gene, "stop"]))+
    theme(axis.text.x = element_text(angle = 90), legend.position = "bottom",
          axis.title.x = element_text(hjust = 0.5, size = 14),
          axis.title.y = element_text(hjust = 0.5, size = 14)
          )
  
  # save plot
  file = paste0(plot_dir, "amplicon_geno_", gene, files_mod)
  svglite(filename = paste0(file, ".svg"), width = 5, height = 5 +(round(sqrt(length(unique(data$sample))))/3*2))
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
    .groups = "drop"
  )

head(stats_summary)

write_tsv(stats_summary, file = paste0(plot_dir, files_mod, "-depth_stats.tsv"))
