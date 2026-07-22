# plot relatedness values from VCFtools --relatedness(2)

# Load Libraries
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  tidyr,
  dplyr,
  gplots,
  tidyverse,
  cowplot,
  pheatmap,
  viridis,
  ggridges,
  readr
)

# --- Define Directories ---
wd <- "C:/Users/NJCB/Documents/Purdue/Pmart_wgs/results/04_relatedness2/"
setwd(wd)

wdII <- "C:/Users/NJCB/Documents/Purdue/Pmart_wgs/figures_and_tables/04_relatedness/"
# Create output directory if it doesn't exist
if (!dir.exists(wdII)) {
  dir.create(wdII, recursive = TRUE)
}

# #############################
# --- Data Loading and Initial Processing ---
# #############################

# Read sample metadata
samples <- read.csv("C:/Users/NJCB/Documents/Purdue/Pmart_wgs/sample_data/Pmart_all_samples.csv")

# Read vcftools --relatedness2 output
relate <- read_tsv("./Pmart_filtered_CHROM_het_filt.relatedness2")

# Read early and late individual lists
early_individuals_df <- read.table("~/Purdue/Pmart_wgs/scripts/E.txt")
late_individuals_df <- read.table("~/Purdue/Pmart_wgs/scripts/L.txt")

# Extract unique individual IDs for early and late groups
early_individuals <- unique(early_individuals_df$V1)
late_individuals <- unique(late_individuals_df$V1)

# Define individual to remove
INDIVIDUAL_TO_REMOVE <- "P_18_53_253"

# --- Consolidated Data Cleaning ---
processed_relate <- relate %>%
  # Remove self-comparisons
  filter(INDV1 != INDV2) %>%
  # Remove the specified individual
  filter(INDV1 != INDIVIDUAL_TO_REMOVE, INDV2 != INDIVIDUAL_TO_REMOVE) %>%
  # Filter based on sample metadata
  filter(INDV1 %in% samples$sample_id, INDV2 %in% samples$sample_id) %>%
  # Ensure unique pairwise comparisons (e.g., A-B and not B-A)
  rowwise() %>%
  mutate(
    ordered_indv1 = min(INDV1, INDV2),
    ordered_indv2 = max(INDV1, INDV2)
  ) %>%
  ungroup() %>%
  distinct(ordered_indv1, ordered_indv2, .keep_all = TRUE) %>%
  select(-ordered_indv1, -ordered_indv2)

# Read and process time-to-death data
time_to_death <- read.csv("C:/Users/NJCB/Documents/Purdue/Pmart_wgs/sample_data/all_sample_hour.csv") %>%
  dplyr::rename(INDV = sample_id, time_to_death_hours = hour)

# Merge time-to-death and calculate delta_time_to_death
relate_merged_with_time <- processed_relate %>%
  left_join(time_to_death, by = c("INDV1" = "INDV")) %>%
  dplyr::rename(time_to_death_1 = time_to_death_hours) %>%
  left_join(time_to_death, by = c("INDV2" = "INDV")) %>%
  dplyr::rename(time_to_death_2 = time_to_death_hours) %>%
  drop_na(time_to_death_1, time_to_death_2) %>%
  mutate(delta_time_to_death = abs(time_to_death_1 - time_to_death_2)) %>%
  # Assign early/late group to each individual
  mutate(
    INDV1_group = case_when(
      INDV1 %in% early_individuals ~ "Early",
      INDV1 %in% late_individuals ~ "Late",
      TRUE ~ "Unknown"
    ),
    INDV2_group = case_when(
      INDV2 %in% early_individuals ~ "Early",
      INDV2 %in% late_individuals ~ "Late",
      TRUE ~ "Unknown"
    )
  ) %>%
  # Create comparison type for coloring
  mutate(
    comparison_type = case_when(
      INDV1_group == "Early" & INDV2_group == "Early" ~ "Early vs. Early",
      (INDV1_group == "Early" & INDV2_group == "Late") |
        (INDV1_group == "Late" & INDV2_group == "Early") ~ "Early vs. Late",
      INDV1_group == "Late" & INDV2_group == "Late" ~ "Late vs. Late",
      TRUE ~ "Other"
    )
  ) %>%
  # Set factor levels for comparison_type
  mutate(
    comparison_type = factor(comparison_type,
                             levels = c("Early vs. Early", "Early vs. Late", "Late vs. Late", "Other")
    ),
    # Ensure RELATEDNESS_PHI is numeric
    RELATEDNESS_PHI = as.numeric(RELATEDNESS_PHI)
  )

# #############################
# --- Plot 1: Heatmap of Relatedness Values ---
# #############################

# Adjust late_individuals for ordering
late_individuals_for_ordering <- late_individuals[!late_individuals == INDIVIDUAL_TO_REMOVE]

# Combine individuals in desired order for heatmap
individuals_in_data <- unique(c(processed_relate$INDV1, processed_relate$INDV2))
ordered_individuals <- c(
  early_individuals[early_individuals %in% individuals_in_data],
  late_individuals_for_ordering[late_individuals_for_ordering %in% individuals_in_data]
)

# Initialize relatedness matrix
relatedness_matrix <- matrix(NA,
                             nrow = length(ordered_individuals),
                             ncol = length(ordered_individuals),
                             dimnames = list(ordered_individuals, ordered_individuals)
)

# Fill relatedness matrix
for (i in 1:nrow(processed_relate)) {
  indv1 <- as.character(processed_relate$INDV1[i])
  indv2 <- as.character(processed_relate$INDV2[i])
  phi_value <- processed_relate$RELATEDNESS_PHI[i]
  
  if (indv1 %in% ordered_individuals && indv2 %in% ordered_individuals) {
    relatedness_matrix[indv1, indv2] <- phi_value
    relatedness_matrix[indv2, indv1] <- phi_value
  }
}

# Create row and column annotations
annotation_row <- data.frame(Group = ifelse(rownames(relatedness_matrix) %in% early_individuals, "Early", "Late"))
rownames(annotation_row) <- rownames(relatedness_matrix)

annotation_col <- data.frame(Group = ifelse(colnames(relatedness_matrix) %in% early_individuals, "Early", "Late"))
rownames(annotation_col) <- colnames(relatedness_matrix)

# Plot heatmap
p_heatmap <- pheatmap(relatedness_matrix,
                      color = viridis(n = 256, alpha = 1, begin = 0, end = 0.9, option = "inferno"),
                      show_rownames = FALSE,
                      show_colnames = FALSE,
                      fontsize_row = 5,
                      fontsize_col = 5,
                      na_col = "lightgrey",
                      cluster_rows = FALSE,
                      cluster_cols = FALSE,
                      annotation_row = annotation_row,
                      annotation_col = annotation_col,
                      annotation_colors = list(Group = c(Early = "lightblue", Late = "darkred"))
)

p_heatmap
# Save heatmap
#ggsave(paste0(wdII, "relatedness2_P-heatmap_excluding_double_sample_het_filt_v2.3.svg"), p_heatmap, dpi = 1000, height = 8, width = 10)
#ggsave(paste0(wdII, "relatedness2_P-heatmap_excluding_double_sample_het_file_v2.3.png"), p_heatmap, dpi = 1000, height = 8, width = 10)


######make an additional plot with relatedness groups to put over half


relat_colors <- c("#000000", "#192bc2", "#0197f6", "#89cff0")

p_heatmap2 <- pheatmap(relatedness_matrix,
                      breaks = c(-Inf, 0.0442, 0.0884, 0.177, Inf),
                      color = relat_colors,
                      show_rownames = FALSE,
                      show_colnames = FALSE,
                      fontsize_row = 5,
                      fontsize_col = 5,
                      na_col = "lightgrey",
                      cluster_rows = FALSE,
                      cluster_cols = FALSE,
                      annotation_row = annotation_row,
                      annotation_col = annotation_col,
                      annotation_colors = list(Group = c(Early = "lightblue", Late = "darkred"),
                      legend = FALSE,
                      annotation_legend = FALSE
                      ))

p_heatmap2

#ggsave(paste0(wdII, "relatedness2_P-heatmap2_excluding_double_sample_het_filt_v2.5.svg"), p_heatmap2, dpi = 1000, height = 8, width = 9.95)


#############################
# --- Plot 2: Histogram of Relatedness Values ---
#############################

# Prepare data for histogram (all unique pairs)
histogram_data <- processed_relate %>%
  mutate(Group = "All samples")

# Plot histogram
h3 <- ggplot(histogram_data, aes(x = RELATEDNESS_PHI, fill = Group)) +
  geom_histogram(bins = 150, fill = "#440154FF") +
  labs(x = "Kinship Coefficient", y = "Counts") +
  theme_cowplot() +
  theme(
    strip.text.x = element_text(size = 14, face = "bold"),
    plot.background = element_rect(fill = "#fffffC"),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 14),
    legend.position = "none",
    panel.grid.major = element_line(color = "gray80", linewidth = 0.5)
  ) +
  scale_x_continuous(breaks = c(-0.1, -0.05, 0, 0.05, 0.1, 0.15, 0.2, 0.25))

h3
# Save main histogram
#ggsave(paste0(wdII, "relatedness2_histogram_excluding_double_sample_het_filt_v2.5.svg"), h3, dpi = 1000, height = 8, width = 10)
#ggsave(paste0(wdII, "relatedness2_histogram_excluding_double_sample_het_filt_v2..png"), h3, dpi = 1000, height = 8, width = 10)

# Create data for the inset plot
inset_data <- histogram_data %>%
  filter(RELATEDNESS_PHI > 0.05)

# Create inset plot
h3_inset <- ggplot(inset_data, aes(x = RELATEDNESS_PHI)) +
  geom_histogram(bins = 30, fill = "#440154FF", color = "white", linewidth = 0.2) +
  theme_cowplot(font_size = 10) +
  theme(
    plot.background = element_rect(color = "navyblue", linewidth = 0.5),
    plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 8),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(color = "white", fill = "white"),
    panel.border = element_rect(color = "navyblue"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  scale_x_continuous(breaks = c(0.05, 0.1, 0.15, 0.2, 0.25))

# Combine main plot and inset plot
h3_with_inset <- ggdraw() +
  draw_plot(h3) +
  draw_plot(h3_inset, x = 0.4, y = 0.45, width = 0.58, height = 0.5)

h3_with_inset
# Save combined histogram with inset
#ggsave(paste0(wdII, "relatedness2_histogram_excluding_double_sample_het_filt_inset_v2.4.svg"), h3_with_inset, dpi = 1000, height = 6, width = 8)
#ggsave(paste0(wdII, "relatedness2_histogram_excluding_double_sample_het_filt_inset_v2.4.png"), h3_with_inset, dpi = 1000, height = 6, width = 8)



#############################
# --- Plot 3: Delta Time to Death vs. Pairwise Relatedness (SI plot) ---
#NOT GOING IN MS (09-26-25)
#############################

# Create bins for Relatedness values
relate_merged_with_time <- relate_merged_with_time %>%
  mutate(relatedness_bin = cut(RELATEDNESS_PHI,
                               breaks = c(-Inf, 0.05, 0.15, 0.25, Inf),
                               labels = c("<0.05", "0.05-0.15", "0.15-0.25", ">0.25"),
                               right = TRUE,
                               include.lowest = TRUE
  )) %>%
  # Set factor levels for relatedness_bin
  mutate(relatedness_bin = factor(relatedness_bin,
                                  levels = c("<0.05", "0.05-0.15", "0.15-0.25", ">0.25")
  ))

# Calculate R-squared and N for each bin (for potential text labels)
r_squared_data <- relate_merged_with_time %>%
  group_by(relatedness_bin) %>%
  summarise(
    r_squared = ifelse(n() > 1, summary(lm(RELATEDNESS_PHI ~ delta_time_to_death))$r.squared, NA),
    n_obs = n()
  ) %>%
  ungroup() %>%
  mutate(
    label = paste0(relatedness_bin, ": R\u00B2 = ", format(r_squared, digits = 2), " (N=", n_obs, ")")
  ) %>%
  filter(!is.na(r_squared))

p2 <- ggplot(relate_merged_with_time, aes(x = RELATEDNESS_PHI, y = delta_time_to_death)) +
  geom_jitter(aes(colour = comparison_type, fill = comparison_type), alpha = 0.25, size = 1, height = 0.25) +
  geom_smooth(method = "lm", se = TRUE, aes(group = relatedness_bin),
              color = "darkgray", linetype = "dashed"
  ) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", linetype = "solid", linewidth = 1) +
  labs(
    x = "Kinship Coefficient",
    y = "Delta Time to Death (Hours)",
    title = "Delta Time to Death vs. Pairwise Relatedness",
    color = "Comparison Type"
  ) +
  theme_cowplot() +
  theme(
    axis.text.y = element_text(size = 10),
    strip.text.x = element_text(size = 12, face = "bold"),
    plot.background = element_rect(fill = "#ffffff"),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 10),
    legend.position = "bottom",
    panel.grid.major = element_line(color = "gray80", linewidth = 0.5)
  ) +
  scale_color_manual(values = c(
    "Early vs. Early" = "#000022",
    "Early vs. Late" = "#FF0000",
    "Late vs. Late" = "#FFE100",
    "Other" = "#000000"
  ))

p2

#ggsave(paste0(wdII, "relatedness2_excluding_double_sample_het_filt_delta_time_to_death_binned.svg"), p2, dpi = 1000, height = 8, width = 10)
#ggsave(paste0(wdII, "relatedness2_excluding_double_sample_het_filt_delta_time_to_death_binned.png"), p2, dpi = 1000, height = 8, width = 10)



#############################
# --- Plot 4: Delta Time to Death vs. Relatedness (Jittered Boxplot) ---
#############################

# Create bins for boxplots
min_phi <- floor(min(relate_merged_with_time$RELATEDNESS_PHI) * 20) / 20
max_phi <- ceiling(max(relate_merged_with_time$RELATEDNESS_PHI) * 20) / 20
phi_breaks <- seq(min_phi, max_phi, by = 0.05)

relate_merged_for_boxplot <- relate_merged_with_time %>%
  mutate(
    relatedness_bin_midpoint = cut(RELATEDNESS_PHI,
                                   breaks = phi_breaks,
                                   labels = (phi_breaks[-length(phi_breaks)] + phi_breaks[-1]) / 2,
                                   right = TRUE,
                                   include.lowest = TRUE,
                                   ordered_result = TRUE
    ) %>% as.character() %>% as.numeric(),
    relatedness_group_for_boxplot = cut(RELATEDNESS_PHI,
                                        breaks = phi_breaks,
                                        labels = seq_along(phi_breaks[-1]),
                                        right = TRUE,
                                        include.lowest = TRUE,
                                        ordered_result = TRUE
    ) %>% as.character()
  )

# Calculate overall mean and median kinship coefficient
mean_phi <- mean(relate_merged_for_boxplot$RELATEDNESS_PHI, na.rm = TRUE)
median_phi <- median(relate_merged_for_boxplot$RELATEDNESS_PHI, na.rm = TRUE)

p3 <- ggplot(relate_merged_for_boxplot, aes(x = RELATEDNESS_PHI, y = delta_time_to_death)) +
  geom_jitter(color = "#1E4165", alpha = 0.25, size = 1, height = 0.25) +
  geom_vline(xintercept = mean_phi, linetype = "solid", color = "darkgreen", linewidth = 0.8) +
  #geom_boxplot(aes(x = relatedness_bin_midpoint, group = relatedness_group_for_boxplot),
               #fill = "grey85",
               #alpha = 0.75,
               #outlier.shape = NA,
               #width = 0.02
  #) +
  scale_x_continuous(breaks = c(-0.1, -0.05, 0, 0.05, 0.1, 0.15, 0.2, 0.25)) +
  labs(
    x = "Kinship Coefficient",
    y=expression(Delta[scriptstyle("Time to Death (Hours)")])
  ) +
  theme_cowplot() +
  theme(
    strip.text.x = element_text(size = 14, face = "bold"),
    plot.background = element_rect(fill = "#fffffC"),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 14),
    legend.position = "none",
    panel.grid.major = element_line(color = "gray80", linewidth = 0.5)
  ) 

p3

# Save jittered boxplot
#ggsave(paste0(wdII, "relatedness2_excluding_double_sample_het_filt_delta_time_to_death_jittered_points_v2.4_SI.svg"), p3, dpi = 1000, height = 6, width = 8)
#ggsave(paste0(wdII, "relatedness2_excluding_double_sample_het_filt_delta_time_to_death_jittered_points_v2.4_SI.png"), p3, dpi = 1000, height = 6, width = 8)


# #############################
# --- Plot 5: Delta Time to Death vs. Relatedness Group (KING Bins) ---
# #############################

# Define KING-based cutoffs:
# Unrelated: <0.0442 | Third-Degree: 0.0442-0.0884 | Second-Degree: 0.0884-0.177 | First-Degree: >0.177

plot5_data <- relate_merged_with_time %>%
  mutate(
    relatedness_group_king = cut(RELATEDNESS_PHI,
                                 breaks = c(-Inf, 0.0442, 0.0884, 0.177, Inf),
                                 labels = c("Unrelated", "Third-Degree", "Second-Degree", "First-Degree"),
                                 right = TRUE,
                                 include.lowest = TRUE
    )
  ) %>%
  mutate(relatedness_group_king = factor(relatedness_group_king,
                                         levels = c("Unrelated", "Third-Degree", "Second-Degree", "First-Degree")
  )) %>%
  drop_na(relatedness_group_king, delta_time_to_death)


# Calculate mean, SD, SE, and number of pairs for each relatedness group
plot5_summary <- plot5_data %>%
  group_by(relatedness_group_king) %>%
  summarise(
    n_pairs = n(),
    mean_delta = mean(delta_time_to_death, na.rm = TRUE),
    sd_delta = sd(delta_time_to_death, na.rm = TRUE),
    se_delta = sd_delta / sqrt(n_pairs),
    max_delta = max(delta_time_to_death, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(n_label = paste0("n=", n_pairs))


#Plot box plots
p5 <- ggplot(plot5_data, aes(x = relatedness_group_king)) + 
  geom_boxplot(aes(y = delta_time_to_death),
               width = 0.5, outlier.shape = NA, fill = "lightblue", color = "black") +
  geom_point(data = plot5_summary, aes(y = mean_delta), color = "red", size = 4, shape = 18) +
  geom_errorbar(data = plot5_summary,
                aes(ymin = mean_delta - se_delta, ymax = mean_delta + se_delta),
                width = 0.15, color = "red", linewidth = 0.8) +
  geom_text(data = plot5_summary,
            aes(y = 10 * 1.05,
                label = n_label),
            vjust = 0, size = 5) +
  labs(
    x = "Kinship Coefficient",
    y = expression(Delta[scriptstyle("Time to Death (Hours)")]),
  ) +
  theme_cowplot() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    panel.grid.major.y = element_line(color = "gray80", linewidth = 0.5)
  )

print(p5)



########################
max_y_for_n_label <- max(plot5_summary$mean_delta + plot5_summary$se_delta, na.rm = TRUE)

p6 <- ggplot(plot5_summary, aes(x = relatedness_group_king, y = mean_delta)) +
  geom_bar(stat = "identity", fill = "#70E0CC", color = "black", linewidth = 0.5) +
    geom_errorbar(aes(ymin = mean_delta - se_delta, ymax = mean_delta + se_delta),
                width = 0.5,
                color = "black", 
                linewidth = 1) +
  geom_text(aes(y = max_y_for_n_label * 1.1,
                label = n_label),
            vjust = 0, size = 5, color = "black") +
  labs(
    x = "Kinship Coefficient",
    y = expression("Mean "* Delta[scriptstyle("Time to Death (Hours)")])
  ) +
  theme_cowplot() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 12),
    panel.grid.major.y = element_line(color = "gray80", linewidth = 0.5)
  )

print(p6)

#ggsave(paste0(wdII, "relatedness2_excluding_double_sample_het_filt_delta_time_to_death_mean_SE_v2.4.svg"), p6, dpi = 1000, height = 6, width = 8)
#ggsave(paste0(wdII, "relatedness2_excluding_double_sample_het_filt_delta_time_to_death_mean_SE_v2.4.png"), p6, dpi = 1000, height = 6, width = 8)

###############################
#test for significance between groups

#Kruskal-wallis test
k_test <- kruskal.test(delta_time_to_death ~ relatedness_group_king, data = plot5_data)
print(k_test)

#then post hoc Mann-Whitney U

pairwise_mwu <- pairwise.wilcox.test(
  x = plot5_data$delta_time_to_death, 
  g = plot5_data$relatedness_group_king,
  p.adjust.method = "BH"  # Use "bonferroni" if you want to be more conservative
)

print(pairwise_mwu)

# Define which groups to compare (you can add or remove pairs here)
my_comparisons <- list( 
  c("Unrelated", "First-Degree"), 
  c("Unrelated", "Second-Degree"),
  c("Unrelated", "Third-Degree") 
)


p6_v2 <- ggplot(plot5_summary, aes(x = relatedness_group_king, y = mean_delta)) +
  geom_bar(stat = "identity", fill = "#70E0CC", color = "black", linewidth = 0.5) +
  geom_errorbar(aes(ymin = mean_delta - se_delta, ymax = mean_delta + se_delta),
                width = 0.5, color = "black", linewidth = 1) +
  # --- Add Statistical Comparisons ---
  stat_compare_means(
    data = plot5_data, # Use the raw data for the test
    aes(x = relatedness_group_king, y = delta_time_to_death),
    comparisons = my_comparisons,
    method = "wilcox.test", # This is the Mann-Whitney U test
    label = "p.signif",     # Options: "p.format" for numbers, "p.signif" for stars
    step.increase = 0.1     # Adjusts the spacing between multiple brackets
  ) +
  labs(
    x = "Kinship Coefficient",
    y = expression("Mean "* Delta[scriptstyle("Time to Death (Hours)")])
  ) +
  theme_cowplot() +
  theme(
    axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
    panel.grid.major.y = element_line(color = "gray80", linewidth = 0.5)
  )

print(p6_v2)

# #############################
# --- Plot 7: Combined Figure (2) (p_heatmap, h3_with_inset, p6) ---
# #############################


p7_combined <- cowplot::plot_grid(
  p_heatmap$gtable, 
  h3_with_inset,   
  p6,               
  ncol = 3,
  labels = c("a", "b", "c"),
  label_fontface = "bold",
  label_fontfamily = "Arial",
  rel_widths = c(1, 1, 1) 
)

print(p7_combined)
