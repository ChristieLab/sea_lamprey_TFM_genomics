library(tidyverse)


wd <- "/scratch/bell/nbackens/Pmart_wgs/vcf_investigation/vep/UGT_invesigation/"
setwd(wd)

vcf_file_path <- "./UGT3.vcf"
sample_info_file_path <- "./Pmart_all_samples.csv"

# 3. Read and process the VCF data first to get sample IDs
vcf_data <- read_tsv(vcf_file_path, comment = "##")

# Get the list of sample IDs present in the VCF file.
vcf_sample_ids <- colnames(vcf_data)[10:ncol(vcf_data)]

# Select relevant columns from VCF and create a new ID if the original is missing
vcf_subset <- vcf_data %>%
  select(`#CHROM`, `POS`, `ID`, `REF`, `ALT`, everything()) %>%
  # Check if the ID column contains the placeholder "." and replace it
  mutate(ID = ifelse(ID == ".", paste0(`#CHROM`, "_", `POS`), ID))

# Reshape the data to a long format
vcf_long <- vcf_subset %>%
  pivot_longer(
    cols = all_of(vcf_sample_ids),
    names_to = "SampleID",
    values_to = "Genotype_Field"
  ) %>%
  mutate(GT = str_extract(Genotype_Field, "^[^:]+"))

# 4. Read and filter sample information data
# The file is a CSV, so we use read_csv()
sample_info <- read_csv(sample_info_file_path) %>%
  # Rename columns to match the script's convention
  rename(SampleID = sample_id, Phenotype = E_or_L) %>%
  # Filter to only include samples present in the VCF file
  filter(SampleID %in% vcf_sample_ids)

# 5. Merge VCF data with filtered sample information
merged_data <- vcf_long %>%
  left_join(sample_info, by = "SampleID") %>%
  # Convert genotypes to a more readable format for plotting
  mutate(
    Genotype = case_when(
      GT == "0/0" ~ "Hom_Ref (AA)",
      GT == "0/1" | GT == "1/0" ~ "Het (Aa)",
      GT == "1/1" ~ "Hom_Alt (aa)",
      TRUE ~ "NA"
    )
  )

# 6. Calculate genotype frequencies for each SNP per phenotype
genotype_freqs <- merged_data %>%
  group_by(ID, Phenotype, Genotype) %>%
  summarise(Count = n(), .groups = "drop") %>%
  # Calculate frequency for each genotype within each group
  group_by(ID, Phenotype) %>%
  mutate(Frequency = Count / sum(Count))

# 7. Create and save a bar chart for each SNP
# The output will be saved to a subdirectory "plots"
dir.create("UGT3_plots", showWarnings = FALSE)

unique_snps <- unique(genotype_freqs$ID)

for (snp_id in unique_snps) {
  plot_data <- genotype_freqs %>%
    filter(ID == snp_id)
  
  p <- ggplot(plot_data, aes(x = Phenotype, y = Frequency, fill = Genotype)) +
    geom_bar(stat = "identity", position = "stack", width = 0.6) +
    geom_text(
      aes(label = scales::percent(Frequency, accuracy = 1)),
      position = position_stack(vjust = 0.5),
      size = 5
    ) +
    labs(
      title = paste("Genotype Frequencies for SNP:", snp_id),
      x = "Phenotype Group",
      y = "Genotype Frequency",
      fill = "Genotype"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom",
      plot.background = element_rect("white")
    ) +
    scale_fill_brewer(palette = "Set2")
  
  ggsave(
    filename = paste0("UGT3_plots/", snp_id, "_genotype_frequencies.png"),
    plot = p,
    width = 8,
    height = 6,
    units = "in",
    dpi = 300
  )
}

print("Plots have been successfully generated and saved to the 'plots' directory.")


