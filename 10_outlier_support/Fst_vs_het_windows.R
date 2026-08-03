#plot of Fst (SNP) vs heterozygosity wihtin high FST windows.
library(readr)
library(ggplot2)

###############
#set up
wd <- "/scratch/bell/nbackens/Pmart_wgs/vcf_investigation/support/"
setwd(wd)
list.files()

###############
#Read in high FST windows in bed format
highlight_regions <- read.table("evl_big_windows.bed", header = FALSE, sep = "\t", 
                                col.names = c("CHROM_ORIG", "REGION_START", "REGION_END"))
#Read in FST for every SNP
evl <- read.table("./Pmart_filtered_CHROM_MLOCI0.8_MAF0.05_HWE0.05_QUAL30_MINDP10_MAXDP50_het_filt_stdFst_EvL.weir.fst", header = TRUE, row.names = NULL)

#clean up
evl <- evl[startsWith(evl$CHROM, "NC"), ]

#Ztransform the values
evl <- evl %>%
  mutate(score = (WEIR_AND_COCKERHAM_FST - mean(WEIR_AND_COCKERHAM_FST, na.rm = TRUE)) / sd(WEIR_AND_COCKERHAM_FST, na.rm = TRUE)) %>%
  select(-WEIR_AND_COCKERHAM_FST) %>% 
  rename(WEIR_AND_COCKERHAM_ZFST = score) # Rename the score column to P

###############
#read in heterozygosity data
het <- read.table("./snp_heterozygosity.hwe", header = FALSE, stringsAsFactors = FALSE)
head(het)

# Rename columns
het <- het %>%
  rename(CHROM = V1, POS = V2, HO = V3)

###############
#scatter plots

# Ensure the POS column in both dataframes is numeric for proper filtering
evl$POS <- as.numeric(evl$POS)
het$POS <- as.numeric(het$POS)

###############
# Merge FST and heterozygosity data frames
merged_data <- inner_join(evl, het, by = c("CHROM", "POS"))

xmin <- min(merged_data$HO)
xmax <- max(merged_data$HO)
ymin <- min(merged_data$WEIR_AND_COCKERHAM_ZFST)
ymax <- max(merged_data$WEIR_AND_COCKERHAM_ZFST)

#check
head(merged_data)
cat(paste0("Total SNPs after merging: ", nrow(merged_data), "\n\n"))

###############
# Generate scatter plots for each high FST region

# Loop through each row (region) in the highlight_regions data frame
for (i in 1:nrow(highlight_regions)) {
  # Extract current region's details
  region_chrom <- highlight_regions$CHROM_ORIG[i]
  region_start <- highlight_regions$REGION_START[i]
  region_end <- highlight_regions$REGION_END[i]
  
  # Filter the merged data to include only SNPs within the current region
  region_data <- merged_data %>%
    filter(CHROM == region_chrom,
           POS >= region_start,
           POS <= region_end)
  
  # Check if there is any data for the current region
  if (nrow(region_data) > 0) {
    # Create a descriptive title for the plot
    plot_title <- paste0("ZFST vs. HO in ", region_chrom, ":", region_start, "-", region_end)
    # Create a file name for saving the plot
    file_name <- paste0("Fst_vs_Het_Region_", region_chrom, "_", region_start, "_", region_end, ".png")
    
    # Generate the scatter plot using ggplot2
    p <- ggplot(region_data, aes(x = HO, y = WEIR_AND_COCKERHAM_ZFST)) +
      geom_point(alpha = 0.4, color = "darkblue") + 
      ylim(ymin, ymax) +
      xlim(xmin, xmax) +
      labs(
        title = plot_title,
        x = expression("H"["O"]),
        y = expression("Weir and Cockerham ZF"["ST"])
      ) +
      theme_minimal() + # Use a minimal theme for a clean look
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"), 
        axis.title = element_text(face = "bold"),
        plot.background = element_rect(fill = 'white')
      )
    
    # Save the plot to the working directory
    ggsave(file_name, plot = p, width = 8, height = 6, dpi = 300)
    cat(paste0("Successfully generated plot: ", file_name, "\n"))
  } else {
    cat(paste0("No data found for region: ", region_chrom, ":", region_start, "-", region_end, ". Skipping plot generation.\n"))
  }
}



