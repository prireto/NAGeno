# visualize seqdepth
suppressPackageStartupMessages({
    library(ggplot2)
    library(hrbrthemes)
    library(readr) # import tsv files
    library(Formula)
    library(ggpubr)
    library(viridis) # pretty colors
    library(scales)
    library(svglite) # save plots as svgs
    library(rlang) # for interpreting strings as symbols (enables pasted strings to be var placeholders)
    library(dplyr)
    library(fuzzyjoin)
    library(showtext)
})

# add narrow font for plotting
font_add_google("Roboto Condensed", "roboto_condensed")
showtext_auto()

args <- commandArgs(trailingOnly = TRUE)
print("Input used for depth plotting:")

files = args[1]
# files = "/home/vera/gueseer/Pipelines/NanoporeAmpliconGenotyping/tutorial/analysis/filtered_bam_sr/depth/"
# files = "/home/gueseer/permanent/App/tools/NAGeno/tutorial/analysis/filtered_bam_sr/depth/" # 3s
print(files)

files_specifier = args[2]
# files_specifier = "SQK-RBK114-24_barcode"
print(files_specifier)

files_mod = args[3]
# files_mod = "_q90_Q30_MAPQ50"
print(files_mod)

# load bc annotation
anno_dir = args[4]
# anno_dir = "/home/vera/gueseer/Pipelines/NanoporeAmpliconGenotyping/tutorial/Src/barcode_assignment.tsv"
# anno_dir = "/home/gueseer/permanent/App/tools/NAGeno/tutorial/Src/barcode_assignment.tsv" # 3s
print(anno_dir)
bc_anno = data.frame(read_tsv(file = anno_dir, col_names = c("sample", "BC")))

# Convert the 'samples' column in anno to a factor with correctly ordered levels
bc_anno$sample = factor(bc_anno$sample, levels = bc_anno$sample[order(as.numeric(sub("S", "", bc_anno$sample)))])


bed_dir = args[5]
# bed_dir = "/home/vera/gueseer/Pipelines/NanoporeAmpliconGenotyping/tutorial/Src/geno_panel_v4.1.bed"
# bed_dir = "/home/gueseer/permanent/App/tools/NAGeno/tutorial/Src/geno_panel_v4.1.bed" # 3s
print(bed_dir)
bed = data.frame(read_tsv(file = bed_dir, col_names = c("chr", "start", "stop", "gene")))

plot_dir = args[6]
# plot_dir = "/home/vera/gueseer/Pipelines/NanoporeAmpliconGenotyping/tutorial/analysis/output/"
# plot_dir = "/home/gueseer/permanent/App/tools/NAGeno/tutorial/analysis/output/depth/" # 3s
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


# dynamically size plot based on n(samples)
# Number of unique samples
n_samples <- length(unique(data$sample))

# Number of columns and rows
ncol <- ceiling(sqrt(n_samples))
nrow <- ceiling(n_samples / ncol)

# Subplot size (inches)
subplot_width <- 2
subplot_height <- 1.5

# Basal space buffer
base_buffer_w <- 1
base_buffer_h <- 3.3

# Total plot size
plot_width <- ncol * subplot_width + base_buffer_w
plot_height <- nrow * subplot_height + base_buffer_h
plot_width
plot_height

# plot depth individually for each gene by sample
for (gene in genes) {
  plot_by_sample = ggplot(subset(data, gene == gene), aes(x = position, y = depth, color = sample))+
    geom_line()+
    theme_ipsum(base_family = "roboto_condensed")+
    ggtitle(paste0(gene, " amplicon - per base depth"))+
    scale_color_manual(values = viridis(length(unique(data$sample))))+
    scale_y_continuous(trans='log10', limits = c(1, 100000), labels = comma)+
    facet_wrap(~sample, ncol = round(sqrt(length(unique(data$sample)))))+
    geom_hline(yintercept = 10, color = "red", linetype = "dashed")+
    xlim(c(bed[bed$gene == gene, "start"], bed[bed$gene == gene, "stop"]))+
    theme(axis.text.x = element_text(angle = 90), legend.position = "none",
          axis.title.x = element_text(hjust = 0.5, size = 14, margin = margin(t = 15)),
          axis.title.y = element_text(hjust = 0.5, size = 14, margin = margin(r = 15))
          )
  
  # save plot
  file = paste0(plot_dir, "Depth_", gene, "_by_sample")
  svglite(filename = paste0(file, ".svg"), width = plot_width, height = plot_height)
  print(plot_by_sample)
  dev.off()
}

# plot depth for each gene in one plot

for (gene in genes) {
  # generate plot
  plot_by_gene = ggplot(subset(data, gene == gene), aes(x = position, y = depth, color = sample))+
    geom_line()+
    theme_ipsum(base_family = "roboto_condensed")+
    ggtitle(paste0(gene, " amplicon - per base depth"))+
    scale_color_manual(values = viridis(length(unique(data$sample))),
                       guide = guide_legend(ncol = pmax(round(sqrt(length(unique(data$sample)))) - round(sqrt(length(unique(data$sample)))/2), 5))
                       )+
    scale_y_continuous(trans='log10', limits = c(1, 100000), labels = comma)+
    geom_hline(yintercept = 10, color = "red", linetype = "dashed")+
    xlim(c(bed[bed$gene == gene, "start"], bed[bed$gene == gene, "stop"]))+
    theme(axis.text.x = element_text(angle = 90), legend.position = "bottom",
          axis.title.x = element_text(hjust = 0.5, size = 14, margin = margin(t = 15)),
          axis.title.y = element_text(hjust = 0.5, size = 14, margin = margin(r = 15))
          )
  
  # save plot
  file = paste0(plot_dir, "Depth_", gene)
  svglite(filename = paste0(file, ".svg"), width = 6, height = 5 +(round(sqrt(length(unique(data$sample))))/3*2))
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

write_tsv(stats_summary, file = paste0(plot_dir, "/Summary_depth_stats.tsv"))
