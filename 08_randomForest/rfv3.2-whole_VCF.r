# Set working directory
setwd("/scratch/bell/nbackens/Pmart_wgs/vcf_investigation/randomForest")

library(genio)
library(randomForest)

# Read chromosome names from file
chroms <- readLines("orig_chroms.txt")

# Read phenotype data
pheno_df <- read.csv("phenotype.csv")

# Read the shared .fam file using the first chromosome as reference
fam_file <- paste0("./Pmart_filtered_CHROM_MLOCI0.8_MAF0.05_HWE0.05_QUAL30_MINDP10_MAXDP50_het_filt_NO_MISSING.fam")
fam <- read_fam(fam_file)

# Match phenotypes to the sample order in the .fam file
pheno_vector <- pheno_df$phenotype[match(fam$id, pheno_df$IID)]
pheno_factor <- factor(pheno_vector, levels = c(0, 1), labels = c("early", "late"))

# Identify samples with non-missing phenotypes
keep <- !is.na(pheno_vector)
pheno_factor <- pheno_factor[keep]


# Read PLINK files
chr_data <- read_plink("./Pmart_filtered_CHROM_MLOCI0.8_MAF0.05_HWE0.05_QUAL30_MINDP10_MAXDP50_het_filt_NO_MISSING")
genotypes <- chr_data$X[keep, ]  # Subset genotypes to samples with phenotypes
genotypes <- t(genotypes)

# Access the SNP information
snp_info <- chr_data$bim

# Access the sample information
sample_info <- chr_data$fam

# --- DEBUGGING STEP: Inspect SNP IDs from .bim file ---
cat("\n--- Debugging SNP IDs from .bim file (snp_info) ---\n")
cat("Column names of snp_info:\n")
print(colnames(snp_info))
cat("\nHead of snp_info (first 5 rows):\n")
print(head(snp_info))
cat("\nAre there any '.' in snp_info$id column (first 10 unique, if any):\n")
print(head(unique(snp_info$id[grep("\\.", snp_info$id)]), 10))
cat("\nNumber of SNPs with '.' as ID:", sum(snp_info$id == "."), "\n")
cat("Total number of SNPs:", nrow(snp_info), "\n")
cat("--- End Debugging SNP IDs ---\n\n")


# --- VERIFICATION STEP: Check dimensions and orientation AFTER TRANSPOSE ---
cat("--- Genotype Matrix Dimensions & Orientation Check (AFTER TRANSPOSE) ---\n")
cat("Dimensions of genotype matrix (rows x columns):", dim(genotypes), "\n")
cat("Expected: Rows = Individuals, Columns = SNPs\n")
cat("Number of individuals (rows):", nrow(genotypes), "\n")
cat("Number of SNPs (columns):", ncol(genotypes), "\n")

# Set row names for individuals (from fam file)
rownames(genotypes) <- sample_info$id

# Set column names for SNPs.
# Use chromosomal position (chr_pos) as SNP identifier if SNP ID is missing/dot
# This is a robust way to create unique identifiers.
if (any(snp_info$id == ".")) {
  warning("Some SNP IDs in .bim file are '.', using CHR:POS as SNP identifier.")
  # Create unique SNP names from Chr and Pos
  snp_names_for_cols <- paste0(snp_info$chr, ":", snp_info$pos)
  # Handle potential duplicate CHR:POS (unlikely for unique SNPs but good practice)
  snp_names_for_cols <- make.names(snp_names_for_cols, unique = TRUE)
  colnames(genotypes) <- snp_names_for_cols
  # Update snp_info to also use these new IDs for consistency in final_results merge
  snp_info$original_id <- snp_info$id # Keep original ID if needed
  snp_info$id <- snp_names_for_cols
} else {
  # If SNP IDs are not '.', use them directly
  colnames(genotypes) <- make.names(snp_info$id, unique = TRUE)
}



# Run randomForest analysis
rf_model <- randomForest(x = genotypes, y = pheno_factor, ntree = 100000, importance = TRUE)

# Extract and sort importance scores
importance_scores <- importance(rf_model)
sorted_importance <- sort(importance_scores[, "MeanDecreaseAccuracy"], decreasing = TRUE)
top_snps <- names(sorted_importance)[1:1000]
top_importance <- sorted_importance[1:1000]

# Create data frame for top SNPs
top_snps_df <- data.frame(SNP = top_snps, Importance = top_importance)

# Sanitize chromosome name for file output (replace '.' with '_')
#chr_name_safe <- gsub("\\.", "_", chr_name)

# Save results to CSV
write.csv(top_snps_df, file = paste0("top_snps_WHOLE_VCF_100k.csv"), row.names = FALSE)

importance_scoresii <- importance(rf_model)
write.csv(importance_scoresii, file = "all_importance_scores_100k.csv", row.names = FALSE)


# Print results to console
cat(sprintf("Chromosome %s Top 10 SNPs:\n", chr_name))
print(top_snps_df)
