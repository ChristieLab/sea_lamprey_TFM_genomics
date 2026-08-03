#Plot Fst, AFD, and RandomForest Importance scores for each chromosome
#Author: NJCB
#Date: 08/31/25
#Data: Weir and Cockerham FST, Allele freq dist. (MC), importance score from RandomForest

library(ggplot2)
library(tidyverse)
library(dplyr)
library(viridis)
library(gridExtra)
library(gtable)
library(svglite)

#####################
##Prepare data
#####################
wd <- "/scratch/bell/nbackens/Pmart_wgs/vcf_investigation/support/"
setwd(wd)
list.files()

## read in AFD 
afd <- read.table("./afds_empirical.txt", quote = "")
colnames(afd) <- c("snp", "chrom", "site", "afd")
afd <- afd %>%
  mutate(across(where(is.character), ~ str_replace_all(.x, "\"", "")))
afd$afd <- as.numeric(afd$afd)
afd$site <- as.numeric(afd$site)

## read in FST
evl <- read.table("./Pmart_filtered_CHROM_MLOCI0.8_MAF0.05_HWE0.05_QUAL30_MINDP10_MAXDP50_het_filt_stdFst_EvL.weir.fst", header=TRUE, row.names = NULL) 

#clean to remove SNPs in unplaced scaffolds
evl <- evl[startsWith(evl$CHROM, "NC"), ]

#Ztransform
evl <- evl %>%
  mutate(score = (WEIR_AND_COCKERHAM_FST - mean(WEIR_AND_COCKERHAM_FST, na.rm = TRUE)) / sd(WEIR_AND_COCKERHAM_FST, na.rm = TRUE))


## read in RF
rf <- read.csv("./all_importance_scores_10k.csv")
#add in chrom and SNP information from the .bim file
bim <- read.table("./Pmart_filtered_CHROM_MLOCI0.8_MAF0.05_HWE0.05_QUAL30_MINDP10_MAXDP50_het_filt_NO_MISSING.bim")
rf <- cbind(rf, bim$V1, bim$V4)

#adjust names
names(rf)[names(rf) == "bim$V1"] <- "chrom"
names(rf)[names(rf) == "bim$V4"] <- "site"
rf <- rf[c("chrom", "site", "early", "late", "MeanDecreaseAccuracy", "MeanDecreaseGini")]

#clean up to save mem
rm(bim)

#write.csv(rf, file = "RandomForest_importance_scores_with_chrom_and_site.csv", quote = FALSE)

## read in outlier_windows
outlier_windows <- read.table("./outlier_windows.txt", header = TRUE)


##Cleanup data to use chrom numbers instead of NCBI accessions
original_values <- c(
  "NC_046069.1", "NC_046070.1", "NC_046071.1", "NC_046072.1", "NC_046073.1",
  "NC_046074.1", "NC_046075.1", "NC_046076.1", "NC_046077.1", "NC_046078.1",
  "NC_046079.1", "NC_046080.1", "NC_046081.1", "NC_046082.1", "NC_046083.1",
  "NC_046084.1", "NC_046085.1", "NC_046086.1", "NC_046087.1", "NC_046088.1",
  "NC_046089.1", "NC_046090.1", "NC_046091.1", "NC_046092.1", "NC_046093.1",
  "NC_046094.1", "NC_046095.1", "NC_046096.1", "NC_046097.1", "NC_046098.1",
  "NC_046099.1", "NC_046100.1", "NC_046101.1", "NC_046102.1", "NC_046103.1",
  "NC_046104.1", "NC_046105.1", "NC_046106.1", "NC_046107.1", "NC_046108.1",
  "NC_046109.1", "NC_046110.1", "NC_046111.1", "NC_046112.1", "NC_046113.1",
  "NC_046114.1", "NC_046115.1", "NC_046116.1", "NC_046117.1", "NC_046118.1",
  "NC_046119.1", "NC_046120.1", "NC_046121.1", "NC_046122.1", "NC_046123.1",
  "NC_046124.1", "NC_046125.1", "NC_046126.1", "NC_046127.1", "NC_046128.1",
  "NC_046129.1", "NC_046130.1", "NC_046131.1", "NC_046132.1", "NC_046133.1",
  "NC_046134.1", "NC_046135.1", "NC_046136.1", "NC_046137.1", "NC_046138.1",
  "NC_046139.1", "NC_046140.1", "NC_046141.1", "NC_046142.1", "NC_046143.1",
  "NC_046144.1", "NC_046145.1", "NC_046146.1", "NC_046147.1", "NC_046148.1",
  "NC_046149.1", "NC_046150.1", "NC_046151.1", "NC_046152.1", "NC_046153.1",
  "NC_001626.1"
)

# Create a lookup table
lookup_table <- data.frame(original_value = original_values, numerical_value = 1:length(original_values))

# Replace the original values with numerical values in data frames
evl$CHROM <- lookup_table$numerical_value[match(evl$CHROM, lookup_table$original_value)]
afd$chrom <- lookup_table$numerical_value[match(afd$chrom, lookup_table$original_value)]
rf$chrom <- lookup_table$numerical_value[match(rf$chrom, lookup_table$original_value)]
outlier_windows$chromosome <- lookup_table$numerical_value[match(outlier_windows$chromosome, lookup_table$original_value)]

#clean for easy viewing
#rm(lookup_table)
#rm(original_values)

#data prepped

#####################
#Prep genome info (i.e. chrom lengths)
#####################
chromosome_lengths <- read.table("Pmar_chrom_size.txt", header = FALSE, col.names = c("chromosome", "length"))
chromosome_lengths <- chromosome_lengths[!startsWith(chromosome_lengths$chromosome, "NW"), ]
chromosome_lengths$chromosome <- lookup_table$numerical_value[match(chromosome_lengths$chromosome, lookup_table$original_value)]

####################
#Plot FST, AFD, and MeanDecreaseGini to each chromosome with non-overlapping windows along the chrom.
#print to pdf.
####################
# Define the window size 
window_size <- 100000

# Calculate the 99th percentile for each metric across ALL data
afd_99th <- quantile(afd$afd, 0.99, na.rm = TRUE)
fst_99th <- quantile(evl$WEIR_AND_COCKERHAM_FST, 0.99, na.rm = TRUE)
rf_99th <- quantile(rf$MeanDecreaseGini, 0.99, na.rm = TRUE) # Using Gini as the importance metric

#Create PDF
pdf("chromosome_plots_100kb_non-overlapping_window.pdf", width = 11, height = 8.5)

#loop through each chromosome
for (chr in chromosome_lengths$chromosome) {
  
  # Filter data for the current chromosome
  # Note: The column names are different, so we must filter each df separately.
  afd_chr <- afd %>% filter(chrom == chr)
  evl_chr <- evl %>% filter(CHROM == chr)
  rf_chr <- rf %>% filter(chrom == chr)
  
  # Bin the data and calculate metrics for each window
  # For afd data
  binned_afd <- afd_chr %>%
    mutate(bin = floor(site / window_size)) %>%
    group_by(bin) %>%
    summarize(
      n_snps = n(),
      mean_afd = mean(afd, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(
      midpoint = (bin * window_size) + (window_size / 2)
    )
  
  # For Fst data
  binned_evl <- evl_chr %>%
    mutate(bin = floor(POS / window_size)) %>%
    group_by(bin) %>%
    summarize(
      n_snps = n(),
      mean_fst = mean(WEIR_AND_COCKERHAM_FST, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(
      midpoint = (bin * window_size) + (window_size / 2)
    )
  
  # For Random Forest importance data
  binned_rf <- rf_chr %>%
    mutate(bin = floor(site / window_size)) %>%
    group_by(bin) %>%
    summarize(
      n_snps = n(),
      mean_rf = mean(MeanDecreaseGini, na.rm = TRUE), # Using Gini as the metric
      .groups = 'drop'
    ) %>%
    mutate(
      midpoint = (bin * window_size) + (window_size / 2)
    )
  
  # Get the total length of the current chromosome
  chr_length <- chromosome_lengths$length[chromosome_lengths$chromosome == chr]
  
  # Create the plot
  p <- ggplot() +
    # Add AFD line and points
    geom_line(data = binned_afd, aes(x = midpoint, y = mean_afd, color = "AFD"), alpha = 0.2) +
    geom_point(data = binned_afd, aes(x = midpoint, y = mean_afd, size = n_snps, color = "AFD")) +
    # Add Fst line and points
    geom_line(data = binned_evl, aes(x = midpoint, y = mean_fst, color = "Fst"), alpha = 0.2) +
    geom_point(data = binned_evl, aes(x = midpoint, y = mean_fst, size = n_snps, color = "Fst")) +
    # Add RF importance line and points
    geom_line(data = binned_rf, aes(x = midpoint, y = mean_rf, color = "Random Forest Importance"), alpha = 0.2) +
    geom_point(data = binned_rf, aes(x = midpoint, y = mean_rf, size = n_snps, color = "Random Forest Importance")) +
    # Add the 99th percentile dotted lines
    geom_hline(yintercept = afd_99th, linetype = "dotted", color = "darkblue", size = 1) +
    geom_hline(yintercept = fst_99th, linetype = "dotted", color = "darkred", size = 1) +
    geom_hline(yintercept = rf_99th, linetype = "dotted", color = "darkgreen", size = 1) +
    
    # Customize plot aesthetics
    scale_size_continuous(range = c(0.5, 5)) + # Adjust point size range
    scale_color_manual(values = c("AFD" = "#440154FF", "Fst" = "#21908CFF", "Random Forest Importance" = "#FDE725FF")) + # Viridis colors
    labs(
      title = paste("Chromosome", chr),
      x = "Genomic Position (bp)",
      y = "Value",
      size = "SNPs in Window",
      color = "Metric"
    ) +
    xlim(0, chr_length) + # Set x-axis limits to the full chromosome length
    theme_classic() +
    theme(legend.position = "bottom")
  
  # Print the plot to the PDF file
  print(p)
}

dev.off()


####################
#Create a multi-panel plot for each chromosome
#using Z-transformed FST values this time.
#Also calculating the 99th percentile of the windowed values rather than all the data
####################
# Define the window size 
window_size <- 100000

# Create an empty list 
binned_data_list <- list()

for (chr in chromosome_lengths$chromosome) {
  # Filter data for the current chromosome
  afd_chr <- afd %>% filter(chrom == chr)
  evl_chr <- evl %>% filter(CHROM == chr)
  rf_chr <- rf %>% filter(chrom == chr)
  
  # Bin the data and calculate metrics for each window
  # For AFD data
  binned_afd <- afd_chr %>%
    mutate(bin = floor(site / window_size)) %>%
    group_by(bin) %>%
    summarize(
      n_snps = n(),
      value = mean(afd, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(midpoint = (bin * window_size) + (window_size / 2), metric = "AFD", chromosome = chr)
  
  # For Fst data
  binned_evl <- evl_chr %>%
    mutate(bin = floor(POS / window_size)) %>%
    group_by(bin) %>%
    summarize(
      n_snps = n(),
      value = mean(score, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(midpoint = (bin * window_size) + (window_size / 2), metric = "Fst (Z-score)", chromosome = chr)
  
  # For Random Forest importance data
  binned_rf <- rf_chr %>%
    mutate(bin = floor(site / window_size)) %>%
    group_by(bin) %>%
    summarize(
      n_snps = n(),
      value = mean(MeanDecreaseGini, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(midpoint = (bin * window_size) + (window_size / 2), metric = "Random Forest Importance", chromosome = chr)
  
  # Combine the three binned data frames for the current chromosome
  combined_chr_data <- bind_rows(binned_afd, binned_evl, binned_rf)
  
  # Add the combined data for this chromosome to the list
  binned_data_list[[as.character(chr)]] <- combined_chr_data
}

# Combine all binned data into a single data frame
all_binned_data <- bind_rows(binned_data_list)

# Calculate the 99th percentile for each metric's binned values
percentiles_99 <- all_binned_data %>%
  group_by(metric) %>%
  summarize(yint = quantile(value, 0.99, na.rm = TRUE), .groups = 'drop')

# Add a 'highlight' column to flag points above the 99th percentile
all_binned_data_highlighted <- all_binned_data %>%
  left_join(percentiles_99, by = "metric") %>%
  mutate(highlight = value > yint)

# Create a new PDF file for these plots
pdf("chromosome_plots_100kb_multi-panel_binned99th.pdf", width = 11, height = 8.5)

# Loop through each chromosome to create the plots using the new data
for (chr in chromosome_lengths$chromosome) {
  # Filter the combined data for the current chromosome
  plot_data <- all_binned_data_highlighted %>% filter(chromosome == chr)
  
  # Get the total length of the current chromosome
  chr_length <- chromosome_lengths$length[chromosome_lengths$chromosome == chr]
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = midpoint, y = value)) +
    # Add a line and point layer for all data
    geom_line(aes(color = metric), alpha = 0.2) +
    geom_point(aes(size = n_snps, color = metric), alpha = 0.5) +
    
    # Add a point layer specifically for the highlighted points
    # The filter is used to ensure we only plot the highlighted points here
    geom_point(data = filter(plot_data, highlight == TRUE), 
               aes(x = midpoint, y = value),
               color = "black", 
               shape = 21,
               size = 4, # Make the highlight point slightly larger
               stroke = 1) + # Add a black outline for contrast
    
    # Add the 99th percentile lines
    geom_hline(data = percentiles_99, aes(yintercept = yint, color = metric), 
               linetype = "dotted", size = 1) +
    
    # Use facet_wrap to create a panel for each metric
    facet_wrap(~ metric, scales = "free_y", ncol = 1) +
    
    # Customize plot aesthetics
    scale_size_continuous(range = c(0.5, 3)) +
    scale_color_manual(values = c("AFD" = "#440154FF", "Fst (Z-score)" = "#21908CFF", "Random Forest Importance" = "#FDE725FF")) +
    labs(
      title = paste("Chromosome", chr),
      x = "Genomic Position (bp)",
      y = "Value",
      size = "SNPs in Window",
      color = "Metric"
    ) +
    xlim(0, chr_length) +
    theme_classic() +
    theme(legend.position = "bottom")
  
  # Print the plot to the PDF file
  print(p)
}

# Close the PDF device
dev.off()


####################
#Look for bins where all metrics are above 95th (conservative) and 99th(strict) percentile
#Make dataframe and plot
####################
# Define the window size for this analysis. You can change this value.
window_size <- 100000

# Create an empty list to store the binned data for each chromosome
binned_data_list <- list()

# Loop through each chromosome and re-bin the data
for (chr in chromosome_lengths$chromosome) {
  # Filter data for the current chromosome
  afd_chr <- afd %>% filter(chrom == chr)
  evl_chr <- evl %>% filter(CHROM == chr)
  rf_chr <- rf %>% filter(chrom == chr)
  
  # Bin the data and calculate metrics for each window
  # For AFD data
  binned_afd <- afd_chr %>%
    mutate(bin = floor(site / window_size)) %>%
    group_by(bin) %>%
    summarize(
      n_snps = n(),
      value = mean(afd, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(midpoint = (bin * window_size) + (window_size / 2), metric = "AFD", chromosome = chr)
  
  # For Fst data
  binned_evl <- evl_chr %>%
    mutate(bin = floor(POS / window_size)) %>%
    group_by(bin) %>%
    summarize(
      n_snps = n(),
      value = mean(score, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(midpoint = (bin * window_size) + (window_size / 2), metric = "Fst (Z-score)", chromosome = chr)
  
  # For Random Forest importance data
  binned_rf <- rf_chr %>%
    mutate(bin = floor(site / window_size)) %>%
    group_by(bin) %>%
    summarize(
      n_snps = n(),
      value = mean(MeanDecreaseGini, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(midpoint = (bin * window_size) + (window_size / 2), metric = "Random Forest Importance", chromosome = chr)
  
  # Combine the three binned data frames for the current chromosome
  combined_chr_data <- bind_rows(binned_afd, binned_evl, binned_rf)
  
  # Add the combined data for this chromosome to the list
  binned_data_list[[as.character(chr)]] <- combined_chr_data
}

# Combine all binned data into a single data frame
all_binned_data <- bind_rows(binned_data_list)

# Calculate the 99th percentile for each metric's binned values
percentiles_99 <- all_binned_data %>%
  group_by(metric) %>%
  summarize(yint = quantile(value, 0.99, na.rm = TRUE), .groups = 'drop')

# Add a 'highlight_99' column to flag points above the 99th percentile
all_binned_data_highlighted <- all_binned_data %>%
  left_join(percentiles_99, by = "metric") %>%
  mutate(highlight_99 = value > yint) %>%
  select(-yint)

# Calculate the 95th percentile for each metric's binned values
percentiles_95 <- all_binned_data %>%
  group_by(metric) %>%
  summarize(yint = quantile(value, 0.95, na.rm = TRUE), .groups = 'drop')

# Add a 'highlight_95' column to flag points above the 95th percentile
all_binned_data_highlighted <- all_binned_data_highlighted %>%
  left_join(percentiles_95, by = "metric") %>%
  mutate(highlight_95 = value > yint) %>%
  select(-yint)

# Identify windows where ALL metrics are above the 99th percentile
top_1_percent_windows <- all_binned_data_highlighted %>%
  group_by(chromosome, bin, midpoint) %>%
  summarize(all_metrics_above_99 = all(highlight_99), .groups = 'drop') %>%
  filter(all_metrics_above_99) %>%
  select(-all_metrics_above_99) %>%
  ungroup() %>%
  mutate(window_size_kb = window_size / 1000)

# Identify windows where ALL metrics are above the 95th percentile
top_5_percent_windows <- all_binned_data_highlighted %>%
  group_by(chromosome, bin, midpoint) %>%
  summarize(all_metrics_above_95 = all(highlight_95), .groups = 'drop') %>%
  filter(all_metrics_above_95) %>%
  select(-all_metrics_above_95) %>%
  ungroup() %>%
  mutate(window_size_kb = window_size / 1000)

# Assign the window size to the data frame names
assign(paste0("top_1_percent_windows_", window_size), top_1_percent_windows)
assign(paste0("top_5_percent_windows_", window_size), top_5_percent_windows)

# Print the resulting data frames
print(paste0("Windows where all metrics are in the top 1% (window size: ", window_size, "):"))
print(get(paste0("top_1_percent_windows_", window_size)))

print(paste0("Windows where all metrics are in the top 5% (window size: ", window_size, "):"))
print(get(paste0("top_5_percent_windows_", window_size)))

##Write data

# Create the filename with the window size for the top 1% data
filename_1 <- paste0("top_1_percent_windows_", window_size / 1000, "kb.csv")

# Write the data frame to a CSV file
write.csv(top_1_percent_windows, filename_1, row.names = FALSE)

# Create the filename with the window size for the top 5% data
filename_5 <- paste0("top_5_percent_windows_", window_size / 1000, "kb.csv")

# Write the data frame to a CSV file
write.csv(top_5_percent_windows, filename_5, row.names = FALSE)

####################
#fig3
####################
#First lets make a Manhattan plot of the AFD scores.
# Define the window size 
window_size <- 10000

# Create an empty list to store the binned data for each chromosome
binned_data_list <- list()

# Loop through each chromosome and re-bin the data
for (chr in chromosome_lengths$chromosome) {
  # Filter data for the current chromosome
  afd_chr <- afd %>% filter(chrom == chr)
  evl_chr <- evl %>% filter(CHROM == chr)
  rf_chr <- rf %>% filter(chrom == chr)
  
  # Bin the data and calculate metrics for each window
  binned_afd <- afd_chr %>%
    mutate(bin = floor(site / window_size)) %>%
    group_by(bin) %>%
    summarize(
      n_snps = n(),
      value = mean(afd, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(midpoint = (bin * window_size) + (window_size / 2), metric = "AFD", chromosome = chr)
  
  binned_evl <- evl_chr %>%
    mutate(bin = floor(POS / window_size)) %>%
    group_by(bin) %>%
    summarize(
      n_snps = n(),
      value = mean(score, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(midpoint = (bin * window_size) + (window_size / 2), metric = "Fst (Z-score)", chromosome = chr)
  
  binned_rf <- rf_chr %>%
    mutate(bin = floor(site / window_size)) %>%
    group_by(bin) %>%
    summarize(
      n_snps = n(),
      value = mean(MeanDecreaseGini, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(midpoint = (bin * window_size) + (window_size / 2), metric = "Random Forest Importance", chromosome = chr)
  
  # Combine the three binned data frames for the current chromosome
  combined_chr_data <- bind_rows(binned_afd, binned_evl, binned_rf)
  binned_data_list[[as.character(chr)]] <- combined_chr_data
}

all_binned_data <- bind_rows(binned_data_list)

# Calculate the cumulative genome position for a Manhattan plot
# We do this for all chromosomes first to ensure correct plotting
cumulative_lengths <- chromosome_lengths %>%
  mutate(cumsum = cumsum(as.numeric(length))) %>%
  select(-length) %>%
  mutate(odd_even = case_when(
    (1:n() %% 2) == 1 ~ "odd",
    TRUE ~ "even"
  ))

# Combine with all binned data to create a genome-wide plot
all_binned_data_plot_all_chr <- all_binned_data %>%
  left_join(cumulative_lengths, by = "chromosome") %>%
  mutate(genome_position = midpoint + (cumsum - chromosome_lengths$length[match(chromosome, chromosome_lengths$chromosome)]))

all_binned_data_plot_all_chr <- all_binned_data_plot_all_chr %>%
  mutate(metric = case_when(
    metric == "Random Forest Importance" ~ "Random Forest",
    TRUE ~ as.character(metric)
  ))

# Create a data frame for custom chromosome labels and their positions
custom_labels <- c(1:10, 15, 20, 25, 30, 40, 50, 60, 70, 85)
axis_df <- cumulative_lengths %>%
  filter(chromosome %in% custom_labels) %>%
  mutate(midpoint = cumsum - (chromosome_lengths$length[match(chromosome, chromosome_lengths$chromosome)] / 2))

# Prepare the 'outlier_windows' data for plotting
outlier_windows_plot <- outlier_windows %>%
  left_join(cumulative_lengths, by = "chromosome") %>%
  mutate(
    start_cum = start + (cumsum - chromosome_lengths$length[match(chromosome, chromosome_lengths$chromosome)]),
    end_cum = end + (cumsum - chromosome_lengths$length[match(chromosome, chromosome_lengths$chromosome)])
  )

# define colors for plots
manhattan_colors <- c("#71706E", "#230F44")
facet_colors <- c("AFD" = "#230F44", "Fst (Z-score)" = "#21908CFF", "Random Forest" = "#47C16EFF")

# --- Plot 1: Manhattan Plot (Genome-Wide AFD) ---

# This plot displays the binned AFD scores across the entire genome.
# Chromosomes are differentiated by color (odd_even), and the X-axis 
# uses calculated cumulative positions for a Manhattan-style visualization.
p1_afd_manhattan <- all_binned_data_plot_all_chr %>%
  filter(metric == "AFD") %>%
  ggplot(aes(x = genome_position, y = value)) +
  geom_point(aes(color = odd_even), size = 0.5, alpha = 0.8) +
  # Highlight defined outlier regions across the genome
  geom_rect(data = outlier_windows_plot, aes(xmin = start_cum, xmax = end_cum, ymin = -Inf, ymax = Inf), fill = "red", alpha = 0.2, inherit.aes = FALSE) +
  scale_color_manual(values = manhattan_colors) +
  scale_x_continuous(
    breaks = axis_df$midpoint,
    labels = axis_df$chromosome
  ) +
  labs(
    x = "Chromosome",
    y = "AFD"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        legend.position = "none")

# Identify all unique chromosomes that contain defined outlier windows.
chromosomes_to_plot <- unique(outlier_windows$chromosome)

metric_labs <- c(
  "AFD" = "AFD",
  "Fst (Z-score)" = expression(paste(italic(F)[italic(ST)], " (Z-score)")),
  "Random Forest" = "Random Forest"
)

metric.labs <- c("AFD", expression(paste(italic(F)[italic(ST)], " (Z-score)")), "Random Forest")
names(metric.labs) <- c("AFD", "Fst (Z-score)", "Random Forest")

for (chr in chromosomes_to_plot) {
  
  current_chr <- as.character(chr)
  
  # Plot 2: Facetted plots for the current chromosome
  p2_faceted <- all_binned_data_plot_all_chr %>%
    filter(chromosome == current_chr) %>% 
    ggplot(aes(x = midpoint, y = value)) +
    geom_rect(data = outlier_windows %>% filter(chromosome == current_chr), 
              aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf), 
              fill = "red", alpha = 0.3, inherit.aes = FALSE) +
    geom_point(aes(color = metric), size = 1.5, alpha = 0.7) +
    geom_line(aes(color = metric), alpha = 0.2) +
    scale_color_manual(
      values = facet_colors,
    ) +
    facet_wrap(~ metric, scales = "free_y", ncol = 1, strip.position = "left",
               labeller = as_labeller(metric.labs)) +
    labs(
      title = paste("Chromosome", current_chr), 
      x = "Genomic Position (bp)",
      y = NULL
    ) +
    theme_classic() +
    theme(legend.position = "none", 
          plot.background = element_rect(fill = "white", colour = NA),
          panel.background = element_rect(fill = "white", colour = NA),
          strip.background = element_blank(),
          strip.placement = "outside",
          strip.text.y = element_text(angle = 90, hjust = 0.5, size = 12),
          axis.title.y = element_blank(),
          axis.text.y = element_text()
    )
  
  # Combine the Manhattan plot (P1) and the chromosome plot (P2)
  # into a single figure using 'grid.arrange'.
  final_layout <- grid.arrange(
    p1_afd_manhattan,
    p2_faceted,
    ncol = 1,
    heights = c(2, 4)
  )
  
  # Construct dynamic file names based on the current chromosome number.
  pdf_filename <- paste0("Fig3_chrom_", current_chr, ".pdf")
  svg_filename <- paste0("Fig3_chrom_", current_chr, ".svg")
  png_filename <- paste0("Fig3_chrom_", current_chr, ".png")
  
  # Save the combined plot to PDF (vector format) and SVG (vector format).
  pdf(pdf_filename, width = 8.5, height = 6.5)
  grid.draw(final_layout)
  dev.off()
  
  svg(svg_filename, width = 8.5, height = 6.5)
  grid.draw(final_layout)
  dev.off()
  
  # Save the combined plot to a PNG
  png(png_filename, width = 8.5, height = 6.5, units = "in", res = 1000)
  grid.draw(final_layout)
  dev.off()
}
