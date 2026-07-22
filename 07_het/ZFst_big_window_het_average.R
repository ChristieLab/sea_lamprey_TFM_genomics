setwd("~/Purdue/Pmart_wgs/results/03_Fst/ZFst_100kb_isolate/het_windows/")

# Read in data
snp_data <- read.table("./snp_heterozygosity.hwe", header = FALSE, stringsAsFactors = FALSE)
window_data <- read.table("evl_big_windows.bed", header = FALSE, stringsAsFactors = FALSE)

# Assign meaningful column names (adjust if your files have headers)
colnames(snp_data) <- c("chr", "pos", "heterozygosity")
colnames(window_data) <- c("chr", "start", "end")

# Initialize a list to store the results
results <- list()
snp_counts <- data.frame(window = character(), n_snps = integer(), stringsAsFactors = FALSE) # Initialize a data frame

# Loop through each window in the window data
for (i in 1:nrow(window_data)) {
  window_chr <- window_data[i, "chr"]
  window_start <- window_data[i, "start"]
  window_end <- window_data[i, "end"]
  window_name <- paste0(window_chr, ":", window_start, "-", window_end) # Create a window identifier
  
  # Subset the SNP data to include only SNPs within the current window
  snps_in_window <- subset(snp_data,
                           chr == window_chr & pos >= window_start & pos <= window_end)
  
  # Calculate the mean heterozygosity and count SNPs
  num_snps <- nrow(snps_in_window)
  
  # Add to the snp_counts data frame
  snp_counts <- rbind(snp_counts, data.frame(window = window_name, n_snps = num_snps))
  
  if (num_snps > 0) {
    mean_heterozygosity <- mean(snps_in_window$heterozygosity)
    results[[i]] <- data.frame(chr = window_chr,
                               start = window_start,
                               end = window_end,
                               mean_heterozygosity = mean_heterozygosity,
                               num_snps = num_snps) 
  } else {
    # If no SNPs in the window
    results[[i]] <- data.frame(chr = window_chr,
                               start = window_start,
                               end = window_end,
                               mean_heterozygosity = NA,
                               num_snps = 0) 
  }
}

# Combine the results into a single data frame
final_results <- do.call(rbind, results)

# Print the results
print(final_results)

# Print the number of SNPs per window
cat("\nNumber of SNPs per window:\n")
print(snp_counts)
