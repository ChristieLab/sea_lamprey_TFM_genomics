#Plot PCA results from PLINK2 for Pmart wgs
#Author: NJCB
#Date: 02/05/2025

#load packages
library(tidyverse)
library(viridis)
library(cowplot)
library(ggsci)
library(ggrepel)
library(grid)
library(gridExtra)
library(Cairo)
library(cowplot)
library(vcfR)
library(ggConvexHull)
library(svglite)
library(readr)
library(dplyr)


#set working directory and check
setwd("C:/Users/NJCB/Documents/Purdue/Pmart_wgs/results/02_PCA/")
getwd()
list.files()
#read in data
pca <- read_table("./Pmart_filtered_CHROM_MLOCI0.8_MAF0.05_HWE0.05_QUAL30_MINDP10_MAXDP50_het_filt_plinkPCA.eigenvec", col_names = TRUE)
eigenval <- scan("./Pmart_filtered_CHROM_MLOCI0.8_MAF0.05_HWE0.05_QUAL30_MINDP10_MAXDP50_het_filt_plinkPCA.eigenval")

#sort out the pca data
#remove nuisance column
pca = select(pca,-1)

# set names
names(pca)[1] <- "ind"
#names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

# sort out the individual species and pops
#read in sample data
samples <- read.csv("C:/Users/NJCB/Documents/Purdue/Pmart_wgs/sample_data/Pmart_all_samples.csv")

#remove unmatched rows from samples file
matching_rows <- samples$sample_id %in% pca$ind
samples_filtered <- samples[matching_rows, ]
samples <- samples_filtered

#remake data frame
pca <- as_tibble(data.frame(pca, samples$actual_length, samples$actual_mass, samples$E_or_L, samples$hour, samples$bag, samples$experiment))

#convert to percentage variance explained.
#Note: use the number of PCs you calculated in pipeline
pve <- data.frame(PC = 1:222, pve = eigenval/sum(eigenval)*100)

#plot percent variance for each PC in a skree
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity", color="black", fill = "dodgerblue2")
a + ylab("Percentage variance explained") + theme_light()
wdII <- ("C:/Users/NJCB/Documents/Purdue/Pmart_wgs/figures_and_tables/02_PCA/")
ggsave(paste(wdII, "Pmart_filtered_CHROM_MLOCI0.8_MAF0.05_HWE0.05_QUAL30_MINDP10_MAXDP50_het_filt_plinkPCA_skree_plot.png", sep = ""), plot = a, height = 5, width = 5, dpi = 1000, units = "in")


#calculate cumulative sum of the percentage variance explained
cumsum(pve$pve)


p1_2 <- ggplot(pca, aes(PC1, PC2, color = samples.E_or_L, shape = samples.experiment)) + 
  geom_point(size = 2, alpha = 0.6) +
  scale_shape_manual(labels = c("1", "2", "3", "4"), values = c(15, 16, 17, 18)) +
  labs(color = "Samples", shape = "Experiment") +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
  scale_color_manual(labels = c("Early", "Late"), values = c("cyan4", "darkred")) +
  theme_cowplot()+
  theme(plot.background = element_rect(fill = "white", color = "white"),
        panel.grid.major.x = element_line(color = "gray", size = 0.5),
        #panel.grid.minor.x = element_line(color = "lightgray", size = 0.25),
        panel.grid.major.y = element_line(color = "gray", size = 0.5))
        #panel.grid.minor.y = element_line(color = "lightgray", size = 0.25))

p1_2


p1_3 <- ggplot(pca, aes(PC1, PC3, color = samples.E_or_L, shape = samples.experiment)) + 
  geom_point(size = 2, alpha = 0.6) +
  scale_shape_manual(labels = c("1", "2", "3", "4"), values = c(15, 16, 17, 18)) +
  labs(color = "Samples", shape = "Experiment") +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC3 (", signif(pve$pve[3], 3), "%)")) +
  theme(panel.background = element_rect(fill = 'white', color = 'white')) +
  scale_color_manual(labels = c("Early", "Late"), values = c("cyan4", "darkred")) +
  theme_cowplot()+
  theme(plot.background = element_rect(fill = "white", color = "white"),
        panel.grid.major.x = element_line(color = "gray", size = 0.5),
        #panel.grid.minor.x = element_line(color = "lightgray", size = 0.25),
        panel.grid.major.y = element_line(color = "gray", size = 0.5))
        #panel.grid.minor.y = element_line(color = "lightgray", size = 0.25))
p1_3

p2_3 <- ggplot(pca, aes(PC2, PC3, color = samples.E_or_L, shape = samples.experiment)) + 
  geom_point(size = 2, alpha = 0.6) +
  scale_shape_manual(labels = c("1", "2", "3", "4"), values = c(15, 16, 17, 18)) +
  labs(color = "Samples", shape = "Experiment") +
  xlab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) + ylab(paste0("PC3 (", signif(pve$pve[3], 3), "%)")) +
  theme(panel.background = element_rect(fill = 'white', color = 'white')) +
  scale_color_manual(labels = c("Early", "Late"), values = c("cyan4", "darkred")) +
  theme_cowplot()+
  theme(plot.background = element_rect(fill = "white", color = "white"),
        panel.grid.major.x = element_line(color = "gray", size = 0.5),
        #panel.grid.minor.x = element_line(color = "lightgray", size = 0.25),
        panel.grid.major.y = element_line(color = "gray", size = 0.5))
#panel.grid.minor.y = element_line(color = "lightgray", size = 0.25))

p2_3


wdII <- ("C:/Users/NJCB/Documents/Purdue/Pmart_wgs/figures_and_tables/02_PCA/")
ggsave(paste(wdII, "Pmart_filtered_CHROM_MLOCI0.8_MAF0.05_HWE0.05_QUAL30_MINDP10_MAXDP50_het_filt_plinkPCA_PC1_PC2.png", sep=""), plot=p1_2, units = "in", dpi=1000, height=6, width=10)
ggsave(paste(wdII, "Pmart_filtered_CHROM_MLOCI0.8_MAF0.05_HWE0.05_QUAL30_MINDP10_MAXDP50_het_filt_plinkPCA_PC1_PC2.svg", sep=""), plot=p1_2, units = "in", dpi=1000, height=6, width=10)
ggsave(paste(wdII, "Pmart_filtered_CHROM_MLOCI0.8_MAF0.05_HWE0.05_QUAL30_MINDP10_MAXDP50_het_filt_plinkPCA_PC1_PC3.png", sep=""), plot=p1_3, units = "in", dpi=1000, height=6, width=10)
ggsave(paste(wdII, "Pmart_filtered_CHROM_MLOCI0.8_MAF0.05_HWE0.05_QUAL30_MINDP10_MAXDP50_het_filt_plinkPCA_PC1_PC3.svg", sep=""), plot=p1_3, units = "in", dpi=1000, height=6, width=10)
ggsave(paste(wdII, "Pmart_filtered_CHROM_MLOCI0.8_MAF0.05_HWE0.05_QUAL30_MINDP10_MAXDP50_het_filt_plinkPCA_PC2_PC3.png", sep=""), plot=p2_3, units = "in", dpi=1000, height=6, width=10)
ggsave(paste(wdII, "Pmart_filtered_CHROM_MLOCI0.8_MAF0.05_HWE0.05_QUAL30_MINDP10_MAXDP50_het_filt_plinkPCA_PC2_PC3.svg", sep=""), plot=p2_3, units = "in", dpi=1000, height=6, width=10)


#####Add in labels to see if outliers have missing data
library(ggrepel)
p1_2 <- ggplot(pca, aes(PC1, PC2, color = samples.E_or_L, shape = samples.experiment, label = ind)) + 
  geom_point(size = 1) +
  geom_text_repel(max.overlaps = 250) +
  scale_shape_manual(labels = c("1", "2", "3", "4"), values = c(15, 16, 17, 18)) +
  labs(color = "Samples", shape = "Experiment") +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
  scale_color_manual(labels = c("Early", "Late"), values = c("cyan4", "darkred")) +
  theme_cowplot()+
  theme(plot.background = element_rect(fill = "white", color = "white"),
        )

p1_2


p1_3 <- ggplot(pca, aes(PC1, PC3, color = samples.E_or_L, shape = samples.experiment, label = ind)) + 
  geom_point(size = 1) +
  geom_text_repel(max.overlaps = 250) +
  scale_shape_manual(labels = c("1", "2", "3", "4"), values = c(15, 16, 17, 18)) +
  labs(color = "Samples", shape = "Experiment") +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC3 (", signif(pve$pve[3], 3), "%)")) +
  theme(panel.background = element_rect(fill = 'white', color = 'white')) +
  scale_color_manual(labels = c("Early", "Late"), values = c("cyan4", "darkred")) +
  theme_cowplot()+
  theme(plot.background = element_rect(fill = "white", color = "white"))

p1_3

p2_3 <- ggplot(pca, aes(PC2, PC3, color = samples.E_or_L, shape = samples.experiment, label = ind)) + 
  geom_point(size = 2) +
  geom_text_repel(max.overlaps = 250) +
  scale_shape_manual(labels = c("1", "2", "3", "4"), values = c(15, 16, 17, 18)) +
  labs(color = "Samples", shape = "Experiment") +
  xlab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) + ylab(paste0("PC3 (", signif(pve$pve[3], 3), "%)")) +
  theme(panel.background = element_rect(fill = 'white', color = 'white')) +
  scale_color_manual(labels = c("Early", "Late"), values = c("cyan4", "darkred")) +
  theme_cowplot()+
  theme(plot.background = element_rect(fill = "white", color = "white"))
p2_3



wdIII <- ("C:/Users/NJCB/Documents/Purdue/Pmart_wgs/figures_and_tables/02_PCA/repel/")
ggsave(paste(wdIII, "Pmart_filtered_CHROM_MLOCI0.8_MAF0.05_HWE0.05_QUAL30_MINDP10_MAXDP50_plinkPCA_PC1_PC2_repel.png", sep=""), plot=p1_2, units = "in", dpi=1000, height=12, width=18)
ggsave(paste(wdIII, "Pmart_filtered_CHROM_MLOCI0.8_MAF0.05_HWE0.05_QUAL30_MINDP10_MAXDP50_plinkPCA_PC1_PC3repel.png", sep=""), plot=p1_3, units = "in", dpi=1000, height=12, width=18)
ggsave(paste(wdIII, "Pmart_filtered_CHROM_MLOCI0.8_MAF0.05_HWE0.05_QUAL30_MINDP10_MAXDP50_plinkPCA_PC2_PC3repel.png", sep=""), plot=p2_3, units = "in", dpi=1000, height=12, width=18)



##################################
########### BACKBURNER ###########
##################################
install.packages("tidyplots")
library(tidyplots)


pca |>
  tidyplot(x = PC1, y = PC2) |> 
  add_data_points()




#Color PCA points by missing data after filtering (lab meeting ?)

#add the data
miss <- read_tsv("./missing_indv_vcftools.imiss")

pca_merged <- left_join(pca, miss[, c("INDV", "F_MISS")], by = c("ind" = "INDV"))

p1_2_miss <- ggplot(pca_merged, aes(PC1, PC2, color = F_MISS, shape = samples.experiment)) + 
  geom_point(size = 2, alpha = 0.75) +
  scale_shape_manual(labels = c("1", "2", "3", "4"), values = c(15, 16, 17, 18)) +
  labs(color = "F_MISS", shape = "Experiment") +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
  scale_color_viridis() +
  theme_cowplot()+
  theme(plot.background = element_rect(fill = "white", color = "white"))

p1_2_miss

wdII <- ("C:/Users/NJCB/Documents/Purdue/Pmart_wgs/figures_and_tables/02_PCA/")
ggsave(paste(wdII, "Pmart_filtered_CHROM_MLOCI0.8_MAF0.05_HWE0.05_QUAL30_MINDP10_MAXDP50_plinkPCA_PC1_PC2_COLORED_BY_IND_MISS.png", sep=""), plot=p1_2_miss, units = "in", dpi=1000, height=6, width=10)
