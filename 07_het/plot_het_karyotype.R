#Plot heterozygosity heatmap with the Fst windows.

library(ggplot2)
library(dplyr)

setwd("~/Purdue/Pmart_wgs/results/07_het/")

# Read in heterozygosity data
snp_data <- read.table("./snp_heterozygosity.hwe", header = FALSE, stringsAsFactors = FALSE)
head(snp_data)

# Rename columns
snp_data <- snp_data %>%
  rename(CHR = V1, SNP = V2, HO = V3)

#modify chromosome names
snp_data <- snp_data %>%
  mutate(CHR = gsub("\\.1$", "", as.character(CHR)))

################################################################
#Read in the chromosome information for the lamprey genome
chromosome_info <- read.csv("~/Purdue/Pmart_wgs/reference/Pmar_chrom_size.csv", header = TRUE)

# Clean chromosome names
chromosome_info$chr <- gsub("\\.1$", "", chromosome_info$chr)

# Create chromosome lengths data frame
chromosome_lengths_df <- chromosome_info %>%
  rename(seqnames = chr, max_bp = length)

print(head(chromosome_lengths_df))

# Find the overall maximum chromosome length
max_length_overall <- max(chromosome_lengths_df$max_bp, na.rm = TRUE)

################################################################
#read in big windows data
evl2 <- read.table("~/Purdue/Pmart_wgs/results/03_Fst/Pmart_filtered_CHROM_MLOCI0.8_MAF0.05_HWE0.05_QUAL30_MINDP10_MAXDP50_het_filt_100kbFst_EvL.windowed.weir.fst", header=TRUE, row.names = NULL)

#clean to remove SNPs in unplaced scaffolds
evl2 <- evl2[startsWith(evl2$CHROM, "NC"), ]

#add in some number
evl2$SNP <- 1:nrow(evl2)

#change column names to work with qqman
colnames(evl2) <- c("CHR", "START", "BP", "N_VARIANTS", "P", "MEAN","SNP")

#zscore
evl2 <- evl2 %>%
  mutate(score = (P - mean(P, na.rm = TRUE)) / sd(P, na.rm = TRUE))

print(evl2)

zfst_threshold <- quantile(evl2$score, 0.99)
zfst_max <- max(evl2$score)
zfst_min <- min(evl2$score)

#Filter for windows above the threshold
evl_top_1_percent <-  evl2[evl2$score >= zfst_threshold, ]

#Filter out chromosomes of interest
#2, 6, 12, 17, 21, 31, 33, 38, 69
keep <- c("NC_046070.1", "NC_046074.1", "NC_046080.1", "NC_046085.1", "NC_046089.1", "NC_046099.1", "NC_046101.1", "NC_046106.1", "NC_046137.1")

evl_top_1_percent_filt <- evl_top_1_percent %>%
  filter(CHR %in% keep)

#Going to find big windows from multiple points above the threhhold
evl_big_windows <- evl_top_1_percent_filt %>%
  group_by(CHR) %>%
  dplyr::summarise(START = min(START),
                   STOP = max(BP)) %>%
  select(CHR, START, STOP) # Select the desired columns in order

#The high values seem to taper off, adding some cushion to capture all overlapping genes.
evl_big_windows <- evl_big_windows %>%
  mutate(START = START - 400000,
         STOP = STOP + 400000)

#prepare big windows for plotting
evl_plot_df <- evl_big_windows %>%
  mutate(seqnames = gsub("\\.1$", "", CHR)) %>%
  left_join(chromosome_lengths_df, by = "seqnames") %>%
  mutate(y_val = match(seqnames, chromosome_lengths_df$seqnames))

################################################################
#Ensure chromosome names in snp_data are consistent
snp_data <- snp_data %>%
  mutate(CHR = gsub("\\.1$", "", as.character(CHR)))

# Create y-positions based on the order in chromosome_lengths_df
chromosomes_all <- chromosome_lengths_df$seqnames
y_positions_ho <- setNames(1:length(chromosomes_all), chromosomes_all)

snp_data <- snp_data %>%
  left_join(chromosome_lengths_df %>% select(seqnames, max_bp), by = c("CHR" = "seqnames")) %>%
  mutate(y_val = y_positions_ho[CHR])

################################################################
ggplot() +
  # Add the chromosome length bars
  geom_rect(data = chromosome_lengths_df,
            aes(xmin = 0, xmax = max_bp, ymin = y_positions_ho[seqnames] - 0.4, ymax = y_positions_ho[seqnames] + 0.4),
            fill = "grey60", alpha = 0.6) +
  geom_segment(data = snp_data,
               aes(x = SNP, xend = SNP, y = y_val - 0.4, yend = y_val + 0.4, color = HO),
               linewidth = 0.2) +
  geom_rect(data = evl_plot_df,
            aes(xmin = START, xmax = STOP, ymin = y_val - 0.3, ymax = y_val + 0.3),
            fill = "chartreuse3", alpha = 0.5) +
  scale_y_continuous(breaks = 1:length(chromosomes_all), labels = 1:length(chromosomes_all)) +
  scale_color_gradient(low = "lightblue", high = "darkblue", name = "Observed\nHeterozygosity") + # Color gradient for HO
  xlab("Genomic Position") +
  ylab("Chromosome") +
  ggtitle("Heterozygosity and Outlier Regions") +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.ticks.y = element_blank()) +
  xlim(0, max_length_overall * 1.05)
