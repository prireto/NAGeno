# visualise seqdepth

# run on de4s

# setwd("/home/vera/gueseer/Projects/cancerna/genotyping/np_amplicon_geno/analysis/geno1_sr-geno-ref/plots")

library(ggplot2)
library(hrbrthemes)
library(readr) # import tsv files
library(ggpubr)
library(viridis) # pretty colors
library(patchwork) # nicer alignment of subplots
library(rmarkdown) 
library(readxl) # import xlsx files
library(ggrepel) # better way of preventing label and text overlaps in ggplot2
library(scales)
library(svglite) # save plots as svgs
library(rlang) # for interpreting strings as symbols ( enables pasted strings to be var placeholders)
library(dplyr)

library(rsvg)


files="/home/vera/gueseer/Projects/cancerna/genotyping/np_amplicon_geno/data/amplicon_geno2-3_filtered_bam_depth/no_trim/sup/"
files_specifier = "SQK-RBK114-24_barcode"
files_mod = "_q90_Q30_MAPQ50"
# files_mod = ""
# for mut cheatsheet
dir = "/home/vera/gueseer/Projects/cancerna/genotyping/analysis/"
plot_dir = "/home/vera/gueseer/Projects/cancerna/genotyping/np_amplicon_geno/data/amplicon_geno2-3_filtered_bam_depth/no_trim/sup/plots/"
man_dir = "/home/vera/gueseer/Projects/cancerna/genotyping/manuscript/depth/no_trim/post_filtering/"
# plot_dir = man_dir
# load bc annotation
bc_anno = data.frame(read_tsv(file = "/home/vera/gueseer/Projects/cancerna/genotyping/np_amplicon_geno/barcode_assignment_geno2-3.tsv", col_names = c("sample", "BC")))

# samples = c("Mel202", "MH081", "MH080", "MH079", "MH037", "MH025", "MH023", "MH022", "MH006", "MH002")

muts = data.frame(read_xlsx(path = paste0(dir, "mut_cheatsheet.xlsx"), col_names = T))
muts = subset(muts, present == "yes")
muts 

# anonymized samples for manuscript
manuscript_ids = data.frame(read_tsv(file = "/home/vera/gueseer/Projects/cancerna/genotyping/manuscript/samplesheet_manuscript.tsv", col_names = T))
colnames(manuscript_ids) = c("sampleID", "sample")
manuscript_ids[manuscript_ids$sample=="x92.1","sample"] = "92.1"
# Convert the 'samples' column to a factor with correctly ordered levels
manuscript_ids$sampleID = factor(manuscript_ids$sampleID, levels = manuscript_ids$sampleID[order(as.numeric(sub("S", "", manuscript_ids$sampleID)))])
head(manuscript_ids)

# table w regions from old mapping startegy w bed file
mut_regions_bed = data.frame(c("SF3B1", "SRSF", "U2AF1_1", "U2AF1_2", "GNA11", "GNAQ"), c(12000, 2200, 1350, 10800, 24000, 78200), c(13500, 3200, 2000, 11800, 25500, 79000), c(2, 17, 21, 21, 19, 9))
colnames(mut_regions_bed) = c("sample", "start", "stop", "chr")
mut_regions_bed

mut_regions_all = data.frame(c("SF3B1", "SRSF2", "GNA11", "GNAQ"), c(197401784, 76736315, 3118362, 77794297), c(197403284, 76737315, 3119862, 77795097), c(2, 17, 19, 9))
colnames(mut_regions_all) = c("sample", "start", "stop", "chr")
mut_regions_all


#########
### Q ###
#########

# use for loop to iterate over all samples defined in bc_anno samplesheet
# maybe this should rathe rbe a list than separate objects, would make all this get adn paste0 crap obsolete
for (samp in 1:length(bc_anno$sample)) {
  # create full file name
  filename = paste0(files_specifier, bc_anno[samp, "BC"], "_", bc_anno[samp, "sample"], files_mod, ".sorted.bam.q.depth.tsv")
  
  print(paste0("Load ", filename, "."))
  
  # load all bam.q.depth.tsv files and create respective df
  # var name is supposed to be sampleQ + add sample col
  assign(paste0(bc_anno[samp, "sample"], "Q"), cbind(data.frame(read_tsv(file = paste0(files, filename), col_names = T)), sample = bc_anno[samp, "sample"]))

  
}


# cat all in one table
dataQ = do.call(rbind, lapply(paste0(bc_anno[, "sample"], "Q"), get))

dataQ$gene="gene"
# old mapping w bed file
# dataQ[dataQ$contig=="chr17:76734115-76737411","gene"]="SRSF2"
# dataQ[dataQ$contig=="chr19:3094362-3123999","gene"]="GNA11"
# dataQ[dataQ$contig=="chr21:43092956-43107578","gene"]="U2AF1"
# dataQ[dataQ$contig=="chr2:197389784-197435093","gene"]="SF3B1"
# dataQ[dataQ$contig=="chr9:77716097-78031811","gene"]="GNAQ"

dataQ[dataQ$contig=="chr17","gene"]="SRSF2"
dataQ[dataQ$contig=="chr19","gene"]="GNA11"
dataQ[dataQ$contig=="chr21","gene"]="U2AF1"
dataQ[dataQ$contig=="chr2","gene"]="SF3B1"
dataQ[dataQ$contig=="chr9","gene"]="GNAQ"

# change x92.1 back to 92.1 - needed to have x to be porperly handled as String before
dataQ[dataQ$sample=="x92.1", "sample"] = "92.1"
unique(dataQ$sample)

# summary(subset(mel202, gene=="GNAQ"))

ggplot(dataQ, aes(x=position, y= depth, fill = sample))+
  geom_line()+
  theme_ipsum()+
  facet_wrap(~gene, scales = "free_x")


#######################
###  GET SEQ RANGES ###
#######################

# get individual np seq ranges per contig and sample => should in general be very consistent but still to be sure

# create seq range df to fill later in for loop - can't be empty bc it looses colnames then, rm first line later
seq_ranges_np = data.frame("contig" = "bla", "gene" = "bla", "start" = -1, "end" = -1, "sample" = "bla")

# create subset of dataQ w only those of min depth 10 for seq range - everything below counts as not sequenced now
dataQ_depth10 = dataQ[dataQ$depth>=10,]

genes = unique(dataQ$gene)
genes
samples = unique(dataQ$sample)
samples
contigs = unique(dataQ$contig)
contigs


# View(dataQ_depth10)
for (sample in samples) {
  for (contig in contigs) {
    # iterate over all positions of this gene, contig and sample to find first and last depth value above 10
    tmp = subset(dataQ_depth10[dataQ_depth10$sample == sample & dataQ_depth10$contig == contig,])
    start = min(tmp$position)
    end = max(tmp$position)
    gene = tmp[tmp$contig == contig, "gene"][1]
    seq_ranges_np = rbind(seq_ranges_np, c(contig, gene, start, end, sample))

  }
  
}
 
head(seq_ranges_np)
seq_ranges_np = seq_ranges_np[-1,]
head(seq_ranges_np)
unique(seq_ranges_np$sample)


#  assign anonymized sample ids
head(dataQ)
head(manuscript_ids)

dataQ = dataQ %>%
  left_join(manuscript_ids[, c("sampleID", "sample")], by = "sample" )

# write_tsv(seq_ranges_np, file = "/home/vera/gueseer/Projects/cancerna/genotyping/np_amplicon_geno/analysis/np_seq_ranges_geno3.tsv")



##########
## MAPQ ##
##########

# # use for loop to iterate over all samples defined in bc_anno samplesheet
# # maybe this should rather be a list than separate objects, would make all this get adn paste0 crap obsolete
# for (samp in 1:length(bc_anno$sample)) {
#   # create full file name
#   filename = paste0(files_specifier, bc_anno[samp, "BC"], "_", bc_anno[samp, "sample"], files_mod, ".sorted.bam.MAPQ.depth.tsv")
#   
#   print(paste0("Load ", filename, "."))
#   
#   # load all bam.q.depth.tsv files and create respective df
#   # var name is supposed to be sampleQ + add sample col
#   assign(paste0(bc_anno[samp, "sample"], "MAPQ"), cbind(data.frame(read_tsv(file = paste0(files, filename), col_names = T)), sample = bc_anno[samp, "sample"]))
#   
#   
# }
# 
# # cat all in one table
# dataMAPQ = do.call(rbind, lapply(paste0(bc_anno[, "sample"], "Q"), get))
# 
# dataMAPQ$gene="gene"
# dataMAPQ[dataMAPQ$contig=="chr17:76734115-76737411","gene"]="SRSF2"
# dataMAPQ[dataMAPQ$contig=="chr19:3094362-3123999","gene"]="GNA11"
# dataMAPQ[dataMAPQ$contig=="chr21:43092956-43107578","gene"]="U2AF1"
# dataMAPQ[dataMAPQ$contig=="chr2:197389784-197435093","gene"]="SF3B1"
# dataMAPQ[dataMAPQ$contig=="chr9:77716097-78031811","gene"]="GNAQ"
# 
# # summary(subset(mel202, gene=="GNAQ"))
# 
# ggplot(dataMAPQ, aes(x=position, y= depth, fill = sample))+
#   geom_line()+
#   theme_ipsum()+
#   facet_wrap(~gene, scales = "free_x")
# 





##############
#### TODO ####
##############
 #  nothing for now



#############
### PLOTS ###
#############

# choose which filter was used => change data and depth_var based on that
filter = "q"
# filter = "MAPQ"



# check depth for each gene by sample
curr_gene = "SF3B1"
sf3b1_by_sample = ggplot(subset(dataQ, gene==curr_gene), aes(x=position, y= depth, color = sampleID))+
  geom_line()+
  # geom_point(data = subset(muts, gene == curr_gene), aes(x = contigPos, y = 0), color = "red", size = 1.5)+
  # geom_text_repel(data = subset(muts, gene == curr_gene), aes(x = contigPos, y = 0, label = mutCDS), hjust = -0.5, color = "red", size = 3, angle = 90)+
  theme_ipsum()+
  ggtitle("SF3B1 amplicon - per base depth")+
  scale_color_manual(values = viridis(length(unique(dataQ$sampleID))))+
  scale_y_continuous(trans='log10', limits = c(1, 100000))+
  facet_wrap(~sampleID)+
  geom_hline(yintercept = 10, color = "red", linetype = "dashed")+
  xlim(c(mut_regions_all[mut_regions_all$sample==curr_gene, "start"], mut_regions_all[mut_regions_all$sample==curr_gene, "stop"]))+
  theme(axis.text.x = element_text(angle = 90), legend.position = "none")
sf3b1_by_sample

file = paste0(plot_dir, "amplicon_geno_", curr_gene, "_by_sample_", files_mod)
svglite(filename = paste0(file, ".svg"), width = 15, height = 10.5)
sf3b1_by_sample
dev.off()
# need this, to make svg editable by inkscape
rsvg_svg(paste0(file, ".svg"), paste0(file, "_plain.svg"))

curr_gene = "SRSF2"
srsf2_by_sample = ggplot(subset(dataQ, gene==curr_gene), aes(x=position, y= depth, color = sampleID))+
  geom_line()+
  # geom_point(data = subset(muts, gene == curr_gene), aes(x = contigPos, y = 1000), color = "red", size = 1.5)+
  # geom_text_repel(data = subset(muts, gene == curr_gene), aes(x = contigPos, y = 1000, label = mutCDS), hjust = -0.5, color = "red", size = 3, angle = 90)+
  theme_ipsum()+
  ggtitle("SRSF2 amplicon - per base depth")+
  scale_color_manual(values = viridis(length(unique(dataQ$sampleID))))+
  scale_y_continuous(trans='log10', limits = c(1, 100000))+
  facet_wrap(~sampleID)+
  geom_hline(yintercept = 10, color = "red", linetype = "dashed")+
  xlim(c(mut_regions_all[mut_regions_all$sample==curr_gene, "start"], mut_regions_all[mut_regions_all$sample==curr_gene, "stop"]))+
  theme(axis.text.x = element_text(angle = 90), legend.position = "none")

file = paste0(plot_dir, "amplicon_geno_", curr_gene, "_by_sample_", files_mod)
svglite(filename = paste0(file, ".svg"), width = 15, height = 10.5)
srsf2_by_sample
dev.off()
# need this, to make svg editable by inkscape
rsvg_svg(paste0(file, ".svg"), paste0(file, "_plain.svg"))

curr_gene = "GNA11"
gna11_by_sample = ggplot(subset(dataQ, gene==curr_gene), aes(x=position, y= depth, color = sampleID))+
  geom_line()+
  # geom_point(data = subset(muts, gene == curr_gene), aes(x = contigPos, y = 1000), color = "red", size = 1.5)+
  # geom_text_repel(data = subset(muts, gene == curr_gene), aes(x = contigPos, y = 1000, label = mutCDS), hjust = -0.5, color = "red", size = 3, angle = 90)+
  theme_ipsum()+
  ggtitle("GNA11 amplicon - per base depth")+
  scale_color_manual(values = viridis(length(unique(dataQ$sampleID))))+
  scale_y_continuous(trans='log10', limits = c(1, 100000))+
  facet_wrap(~sampleID)+
  geom_hline(yintercept = 10, color = "red", linetype = "dashed")+
  xlim(c(mut_regions_all[mut_regions_all$sample==curr_gene, "start"], mut_regions_all[mut_regions_all$sample==curr_gene, "stop"]))+
  theme(axis.text.x = element_text(angle = 90), legend.position = "none")

file = paste0(plot_dir, "amplicon_geno_", curr_gene, "_by_sample_", files_mod)
svglite(filename = paste0(file, ".svg") , width = 15, height = 10.5)
gna11_by_sample
dev.off()
# need this, to make svg editable by inkscape
rsvg_svg(paste0(file, ".svg"), paste0(file, "_plain.svg"))

curr_gene = "GNAQ"
gnaq_by_sample = ggplot(subset(dataQ, gene==curr_gene), aes(x=position, y= depth, color = sampleID))+
  geom_line()+
  # geom_point(data = subset(muts, gene == curr_gene), aes(x = contigPos, y = 1000), color = "red", size = 1.5)+
  # geom_text_repel(data = subset(muts, gene == curr_gene), aes(x = contigPos, y = 1000, label = mutCDS), hjust = -0.5, color = "red", size = 3, angle = 90)+
  theme_ipsum()+
  ggtitle("GNAQ amplicon - per base depth")+
  scale_color_manual(values = viridis(length(unique(dataQ$sampleID))))+
  scale_y_continuous(trans='log10', limits = c(1, 100000))+
  facet_wrap(~sampleID)+
  geom_hline(yintercept = 10, color = "red", linetype = "dashed")+
  xlim(c(mut_regions_all[mut_regions_all$sample==curr_gene, "start"], mut_regions_all[mut_regions_all$sample==curr_gene, "stop"]))+
  theme(axis.text.x = element_text(angle = 90), legend.position = "none")

file = paste0(plot_dir, "amplicon_geno_", curr_gene, "_by_sample_", files_mod)
svglite(filename = paste0(file, ".svg"), width = 15, height = 10.5)
gnaq_by_sample
dev.off()
#  need this to make editable in inkscape
rsvg_svg(paste0(file, ".svg"), paste0(file, "_plain.svg"))

############################################################################
# plot all samples in one

# set data and depth_var based on previously picked filter (q or MAPQ)
if(filter == "q"){
  data = dataQ
} else if (filter == "MAPQ"){
  data = data_MAPQ
}


# for qAll use this and don't run for loop
depth_var = rlang::sym("depth")

# use depth_var to iterate over different depth cols - !! ensures a string is interpreted as var
# create single plots for each gene using the gene var curr_gene
curr_gene = "GNAQ"
gnaq = ggplot(subset(data, gene==curr_gene), aes(x=position, y= !!depth_var, color = sampleID))+
  geom_line()+
  geom_point(data = subset(muts, gene == curr_gene), aes(x = mutLocus, y = 10), color = "red", size = 1.5)+
  geom_text_repel(data = subset(muts, gene == curr_gene), aes(x = mutLocus, y = 10, label = mutCDS), hjust = 1.5, color = "red", size = 3, angle = 90)+
  theme_ipsum()+
  ylim(0,10000)+
  geom_hline(yintercept = 10, color = "red", linetype = "dashed")+
  # ggtitle(paste0("GNAQ amplicon - per base depth"))+
  # ggtitle("GNAQ amplicon - chr 9")+
  scale_color_manual(values = viridis(length(unique(data$sampleID))))+
  scale_y_continuous(trans='log10', limits = c(1, 100000))+
  scale_x_continuous(breaks = scales::pretty_breaks(n=2),
                     limits = c(mut_regions_all[mut_regions_all$sample==curr_gene, "start"], mut_regions_all[mut_regions_all$sample==curr_gene, "stop"]))+
  labs(subtitle = "GNAQ amplicon - chr 9")+
  theme(
    # axis.title.x = element_text(hjust = 0.5, size = 13, margin = margin(t = 12)),
    axis.title.x = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    plot.subtitle = element_text(size = 16),
    axis.title.y = element_text(hjust = 0.5, size = 15, margin = margin(r = 12))
    )+
  guides(color = guide_legend(nrow = 3))
  

file = paste0(plot_dir, "amplicon_geno", files_mod, "_", curr_gene)
# svglite(filename = paste0(file, ".svg"), width = 7, height = 5)
gnaq
# dev.off()
#  need this to make editable in inkscape
# rsvg_svg(paste0(file, ".svg"), paste0(file, "_plain.svg"))

curr_gene = "GNA11"
gna11 = ggplot(subset(data, gene==curr_gene), aes(x=position, y= !!depth_var, color = sampleID))+
  geom_line()+
  geom_point(data = subset(muts, gene == curr_gene), aes(x = mutLocus, y = 10), color = "red", size = 1.5)+
  geom_text_repel(data = subset(muts, gene == curr_gene), aes(x = mutLocus, y = 10, label = mutCDS), hjust = 1.5, color = "red", size = 3, angle = 90)+
  scale_color_manual(values = viridis(length(unique(data$sampleID))))+
  # ggtitle(paste0("GNA11 amplicon - per base depth"))+
  # ggtitle("GNA11 amplicon - chr 19")+
  theme_ipsum()+
  geom_hline(yintercept = 10, color = "red", linetype = "dashed")+
  scale_y_continuous(trans='log10', limits = c(1, 100000))+
  scale_x_continuous(breaks = scales::pretty_breaks(n=3),
                     # labels = waiver(),
                     limits = c(mut_regions_all[mut_regions_all$sample==curr_gene, "start"], mut_regions_all[mut_regions_all$sample==curr_gene, "stop"]))+
  labs(subtitle = "GNA11 amplicon - chr 19")+
  theme(
    # axis.title.x = element_text(hjust = 0.5, size = 13, margin = margin(t = 12)),
      legend.position = "bottom",
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 15),
      plot.subtitle = element_text(size = 16),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
      # axis.title.y = element_text(hjust = 0.5, size = 13, margin = margin(r = 12))
      )+
  guides(color = guide_legend(nrow = 3))

file = paste0(plot_dir, "amplicon_geno", files_mod, "_", curr_gene)
# svglite(filename = paste0(file, ".svg"), width = 7, height = 5)
gna11
# dev.off()
#  need this to make editable in inkscape
# rsvg_svg(paste0(file, ".svg"), paste0(file, "_plain.svg"))

curr_gene = "SF3B1"
sf3b1 = ggplot(subset(data, gene==curr_gene), aes(x=position, y= !!depth_var, color = sampleID))+
  geom_line()+
  geom_point(data = subset(muts, gene == curr_gene), aes(x = mutLocus, y = 10), color = "red", size = 1.5)+
  geom_text_repel(data = subset(muts, gene == curr_gene), aes(x = mutLocus, y = 10, label = mutCDS), hjust = 1.5, color = "red", size = 3, angle = 90)+
  theme_ipsum()+
  geom_hline(yintercept = 10, color = "red", linetype = "dashed")+
  # ggtitle(paste0("SF3B1 amplicon - per base depth"))+
  # ggtitle("SF3B1 amplicon - chr 2")+
  scale_color_manual(values = viridis(length(unique(data$sampleID))))+
  scale_y_continuous(trans='log10', limits = c(1, 100000))+
  scale_x_continuous(limits = c(mut_regions_all[mut_regions_all$sample==curr_gene, "start"], mut_regions_all[mut_regions_all$sample==curr_gene, "stop"]),
                     breaks = scales::pretty_breaks(n=2)
                     # labels = label_number(suffix = "")
                     # labels = scales::label_scientific()
                     )+
  labs(subtitle = "SF3B1 amplicon - chr 2")+
  theme(
    # axis.title.x = element_text(hjust = 0.5, size = 13, margin = margin(t = 12)),
    legend.position = "bottom",
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    plot.subtitle = element_text(size = 16),
    axis.title.y = element_text(hjust = 0.5, size = 15, margin = margin(r = 12))
    )+
  guides(color = guide_legend(nrow = 3))

file = paste0(plot_dir, "amplicon_geno", files_mod, "_", curr_gene)
# svglite(filename = paste0(file, ".svg"), width = 7, height = 5)
sf3b1
# dev.off()
#  need this to make editable in inkscape
# rsvg_svg(paste0(file, ".svg"), paste0(file, "_plain.svg"))

curr_gene = "SRSF2"
srsf2 = ggplot(subset(data, gene==curr_gene), aes(x=position, y= !!depth_var, color = sampleID))+
  geom_line()+
  geom_point(data = subset(muts, gene == curr_gene), aes(x = mutLocus, y = 10), color = "red", size = 1.5)+
  geom_text_repel(data = subset(muts, gene == curr_gene), aes(x = mutLocus, y = 10, label = mutCDS), hjust = 1.5, color = "red", size = 3, angle = 90)+
  scale_color_manual(values = viridis(length(unique(data$sampleID))))+
  # ggtitle(paste0("SRSF2 amplicon - per base depth"))+
  # ggtitle("SRSF2 amplicon - chr 17")+
  theme_ipsum()+
  geom_hline(yintercept = 10, color = "red", linetype = "dashed")+
  scale_y_continuous(trans='log10', limits = c(1, 100000))+
  scale_x_continuous(breaks = scales::pretty_breaks(n=2),
                     limits = c(mut_regions_all[mut_regions_all$sample==curr_gene, "start"], mut_regions_all[mut_regions_all$sample==curr_gene, "stop"]))+
  labs(subtitle=("SRSF2 amplicon - chr 17"))+
  theme(
    # axis.title.x = element_text(hjust = 0.5, size = 13, margin = margin(t = 12)),
    legend.position = "bottom",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.subtitle = element_text(size = 16),
    axis.text.x = element_text(size = 15),
    plot.title = element_text(vjust = 0.5),
    axis.text.y = element_blank()
    # axis.title.y = element_text(hjust = 0.5, size = 13, margin = margin(r = 12))
    )+
  guides(color = guide_legend(nrow = 3))

file = paste0(plot_dir, "amplicon_geno", files_mod, "_", curr_gene)
# svglite(filename = paste0(file, ".svg"), width = 7, height = 5)
srsf2
# dev.off()
#  need this to make editable in inkscape
# rsvg_svg(paste0(file, ".svg"), paste0(file, "_plain.svg"))

# save individual plots w ylim -150, and combined plot w ylim XXX to make label not overlap w geom_lines

# include guie and legend.pos lines in previous plots for combined plot
depths_plot = ggarrange(gnaq, gna11, sf3b1, srsf2, common.legend = T, legend = "bottom", nrow = 2, ncol = 2, widths = c(1.1, 1))
# save as 1000x700
depths_plot



# save plot - also save once as 10x10 for legend
file = paste0(plot_dir, "amplicon_geno", files_mod)
svglite(filename = paste0(file, ".svg"), width = 7.5, height = 7.5)
print(annotate_figure(
  depths_plot,
  top = text_grob("Post filtering", face = "bold", size = 16, vjust = 1)
)) # print ensures the plot is actually printed, otherwise timeout before saving
dev.off()
#  need this to make editable in inkscape
rsvg_svg(paste0(file, ".svg"), paste0(file, "_plain.svg"))



# save plot - also save once as 10x10 for legend
file = paste0(plot_dir, "amplicon_geno", files_mod)
svglite(filename = paste0(file, "_legend.svg"), width = 10, height = 10)
print(depths_plot) # print ensures the plot is actually printed, otherwise timeout before saving
dev.off()
#  need this to make editable in inkscape
rsvg_svg(paste0(file, "_legend.svg"), paste0(file, "_plain_legend.svg"))

#  test
ggarrange(blank, blank) %>%
  annotate_figure(text_grob("Title", face = "bold", size = 16))


###################
### DEPTH STATS ###
###################

# use all depth data from 10bp in to avoid the frizzy ends - ugly hardcoded
# take all positions after first in range with depth>1 and before last in range w depth>1
min_2 = min(subset(dataQ, contig == "chr2" & position > mut_regions_all[mut_regions_all$chr == "2", "start"] & depth > 1)$position) + 0
min_9 = min(subset(dataQ, contig == "chr9" & position > mut_regions_all[mut_regions_all$chr == "9", "start"] & depth > 1)$position) + 0
min_17 = min(subset(dataQ, contig == "chr17" & position > mut_regions_all[mut_regions_all$chr == "17", "start"] & depth > 1)$position) + 0
min_19 = min(subset(dataQ, contig == "chr19" & position > mut_regions_all[mut_regions_all$chr == "19", "start"] & depth > 1)$position) + 0

max_2 = max(subset(dataQ, contig == "chr2" & position < mut_regions_all[mut_regions_all$chr == "2", "stop"] & depth > 1)$position) - 0
max_9 = max(subset(dataQ, contig == "chr9" & position < mut_regions_all[mut_regions_all$chr == "9", "stop"] & depth > 1)$position) - 0
max_17 = max(subset(dataQ, contig == "chr17" & position < mut_regions_all[mut_regions_all$chr == "17", "stop"] & depth > 1)$position) - 0
max_19 = max(subset(dataQ, contig == "chr19" & position < mut_regions_all[mut_regions_all$chr == "19", "stop"] & depth > 1)$position) - 0

max_2 - min_2
max_17 - min_17
max_19 - min_19
max_9 - min_9



data_stats = subset(dataQ, 
                    contig == "chr2" & position > min_2 & position < max_2 |
                      contig == "chr9" & position > min_9 & position < max_9 |
                      contig == "chr17" & position > min_17 & position < max_17 |
                      contig == "chr19" & position > min_19 & position < max_19
                    )

stats_all <- data_stats %>%
  select(sampleID, sample, contig, depth, gene) %>%
  group_by(sampleID, sample, contig, gene) %>%
  summarise(
    median = median(depth, na.rm = TRUE),
    mean = round(mean(depth, na.rm = TRUE), 2),
    # min = min(depth, na.rm = TRUE),
    # max = max(depth, na.rm = TRUE),
    .groups = "drop"
  )

stats_sample <- data_stats %>%
  select(sampleID, sample, contig, depth, gene) %>%
  group_by(sampleID, sample) %>%
  summarise(
    median = median(depth, na.rm = TRUE),
    mean = mean(depth, na.rm = TRUE),
    min = min(depth, na.rm = TRUE),
    max = max(depth, na.rm = TRUE),
    .groups = "drop"
  )

stats_contig <- data_stats %>%
  select(sampleID, sample, contig, depth, gene) %>%
  group_by(contig, gene) %>%
  summarise(
    median = median(depth, na.rm = TRUE),
    mean = mean(depth, na.rm = TRUE),
    min = min(depth, na.rm = TRUE),
    max = max(depth, na.rm = TRUE),
    .groups = "drop"
  )

median_overall = median(data_stats$depth)
median_overall

mean_overall = mean(data_stats$depth)
mean_overall

median_g11 = median(subset(data_stats, gene == "GNA11")$depth)
median_g11

median_gq = median(subset(data_stats, gene == "GNAQ")$depth)
median_gq

median_sf = median(subset(data_stats, gene == "SF3B1")$depth)
median_sf

median_sr = median(subset(data_stats, gene == "SRSF2")$depth)
median_sr



write_tsv(stats_all, file = paste0(plot_dir, mod, "-depth_stats.tsv"))
