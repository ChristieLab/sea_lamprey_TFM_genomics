#---------------- Data Description ------------------------------------------
# script written by allison nalesnik, last revised September 13, 2024
# modified by NJCB 5/15/25
# data consists of samples from lamprey experiments at HBBS in August 2022

# Set working directory and load packages
setwd("~/Purdue/Pmart_wgs/results/09_survival/")
list.files()
library(data.table)
library(dplyr)
library(viridis) # Although viridis is loaded, it's not explicitly used in the provided plots. Keep for potential future use.
library(ggplot2)
library(patchwork)

# Define common save directory and plot version
save_dir <- "~/Purdue/Pmart_wgs/figures_and_tables/09_survival/"
plot_version <- "v2.2"

#---------------------------------------- Data Loading and Initial Processing -----------------------------------------------------
# Load data
all_fish <- read.csv("./all_experiment_fish1.csv", header = TRUE) # all experiment fish
samples <- read.csv("./novogene_official.csv", header = TRUE) # all sequenced fish

# Extract names of sequenced fish
sequenced_samples <- samples$fish_ID # Assuming fish_ID is the correct column based on subsequent usage

# Add column to all_fish designating if sequenced or preserved
all_fish$status <- ifelse(all_fish$fish_ID %in% sequenced_samples, "sequenced", "preserved")

# Filter out experiment 2 tank S4
all_fish <- all_fish %>%
  filter(!(experiment == "experiment_2" & tank == "S4"))

# Add statusII for main text figure (early vs extra_late)
all_fish <- all_fish %>%
  mutate(statusII = ifelse(hour > 12, "extra_late", "average"))
print(unique(all_fish$statusII))

# Prepare data for length distribution plots
filtered_fish <- all_fish %>%
  filter(
    !(experiment == "experiment_1" & tank == "S4" & hour %in% c(4, 5)),
    !(experiment == "experiment_1" & tank == "S6" & hour %in% c(4, 5)),
    !(experiment == "experiment_2" & tank == "S6" & hour %in% c(4, 5)),
    !(experiment == "experiment_2" & tank == "S4" & hour %in% c(4, 5)),
    !(experiment == "experiment_3" & tank == "S4" & hour %in% c(5, 6)),
    !(experiment == "experiment_3" & tank == "S6" & hour %in% c(4, 5, 6)),
    !(experiment == "experiment_4" & tank == "S4" & hour == 5),
    !(experiment == "experiment_4" & tank == "S6" & hour %in% c(5, 6))
  )

# Create 'time_group' column based on 'hour' for length distribution plots
all_fish_modified <- filtered_fish %>%
  mutate(time_group = ifelse(hour <= 5, "Early", "Late"))

# Ensure 'time_group' is a factor with desired order
all_fish_modified$time_group <- factor(all_fish_modified$time_group,
                                       levels = c("Early", "Late"))

# Filter data for only "sequenced" samples for the second length distribution plot
sequenced_fish_modified <- all_fish_modified %>%
  filter(status == "sequenced")

# Define common plot parameters for length distributions
x_min_val_length <- floor(min(all_fish_modified$total_length) / 10) * 10 - 10
x_max_val_length <- ceiling(max(all_fish_modified$total_length) / 10) * 10 + 10
x_breaks_val_length <- seq(40, 140, by = 20)

# Define common color scheme for length distributions
common_colors_length_fill <- c("Early" = "#56B4E9", # Light blue for fill
                               "Late" = "#E69F00") # Orange for fill

common_colors_length_mean <- c("Early" = "darkblue", # Dark blue for mean line
                               "Late" = "darkred")  # Dark red for mean line

#############################
# --- Plot 1: Survival Time vs. Length (Sequenced Highlighted) ---
#############################

# Create plot with sequenced samples highlighted
p <- ggplot(all_fish, aes(x = hour, y = total_length)) +
  geom_jitter(aes(color = status, fill = status),
              width = 0.4, size = 1.5, alpha = 0.6) + # Apply alpha to all jitter points
  geom_boxplot(data = subset(all_fish, hour >= 1 & hour <= 10),
               aes(group = hour),
               fill = NA,
               color = "black",
               width = 0.5,
               outlier.shape = NA) +
  scale_color_manual(values = c("sequenced" = "#1E4165", "preserved" = "gray60")) +
  scale_fill_manual(values = c("sequenced" = "#1E4165", "preserved" = "gray60")) + # Also define fill for consistency
  scale_x_continuous(breaks = 1:18) +
  scale_y_continuous(breaks = c(40, 60, 80, 100, 120, 140)) +
  labs(x = "Survival time (hours in TFM)", y = "Total Length (mm)") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        panel.grid.major.y = element_line(),
        panel.background = element_rect("white"))

print(p)

# Save plot
ggsave(paste0(save_dir, "Pmar_tfm_survival_length_hour_boxplot_sequencedSI_", plot_version, ".png"), p, dpi = 1000, height = 6, width = 8)
ggsave(paste0(save_dir, "Pmar_tfm_survival_length_hour_boxplot_sequencedSI_", plot_version, ".svg"), p, dpi = 1000, height = 6, width = 8)

#############################
# --- Plot 2: Survival Time vs. Length (Extra Late Highlighted) ---
#############################

# Create plot with extra late samples highlighted
p1 <- ggplot(all_fish, aes(x = hour, y = total_length)) +
  geom_jitter(aes(color = statusII, fill = statusII),
              width = 0.4, size = 1.5, alpha = 0.6) + # Alpha applied to all points
  geom_boxplot(data = subset(all_fish, hour >= 1 & hour <= 10),
               aes(group = hour),
               fill = NA,
               color = "black",
               width = 0.5,
               outlier.shape = NA) +
  scale_color_manual(values = c("extra_late" = "#B00000", "average" = "gray60")) +
  scale_fill_manual(values = c("extra_late" = "#B00000", "average" = "gray60")) + # Also define fill for consistency
  scale_x_continuous(breaks = 1:18) +
  scale_y_continuous(breaks = c(40, 60, 80, 100, 120, 140)) +
  labs(x = "Survival time (hours in TFM)", y = "Total Length (mm)") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        panel.grid.major.y = element_line(),
        panel.background = element_rect("white"),
        legend.position = "none")

print(p1)

# Save plot
ggsave(paste0(save_dir, "Pmar_tfm_survival_length_hour_boxplot_sequenced_MAIN_", plot_version, ".png"), p1, dpi = 1000, height = 6, width = 8)
ggsave(paste0(save_dir, "Pmar_tfm_survival_length_hour_boxplot_sequenced_MAIN_", plot_version, ".svg"), p1, dpi = 1000, height = 6, width = 8)


#############################
# --- Plot 3: Size Distributions of Early and Late Samples (Before/After Size Matching) ---
#############################

# Calculate max counts for placing annotations (adjust y-position dynamically)
hist_data_all <- hist(all_fish_modified$total_length, breaks=seq(x_min_val_length, x_max_val_length, by=2), plot = FALSE)
y_anno_all <- max(hist_data_all$counts) * 0.95

hist_data_sequenced <- hist(sequenced_fish_modified$total_length, breaks=seq(x_min_val_length, x_max_val_length, by=2), plot = FALSE)
y_anno_sequenced <- max(hist_data_sequenced$counts) * 0.95

# Calculate mean for each time_group for all samples
summary_stats_all <- all_fish_modified %>%
  group_by(time_group) %>%
  summarise(
    mean_length = mean(total_length, na.rm = TRUE)
  )

# Calculate mean for each time_group for sequenced samples
summary_stats_sequenced <- sequenced_fish_modified %>%
  group_by(time_group) %>%
  summarise(
    mean_length = mean(total_length, na.rm = TRUE)
  )

preserved_fish <- subset(all_fish_modified, status == "preserved")
t.test(total_length ~ time_group, data = preserved_fish)
#p < 2.2e-16
sequenced_fish <- subset(all_fish_modified, status == "sequenced")
t.test(total_length ~ time_group, data = sequenced_fish)
#p = 1.652e-05

# All Samples plot
p_all_counts <- ggplot(all_fish_modified, aes(x = total_length, fill = time_group)) +
  geom_histogram(binwidth = 2,
                 color = "black",
                 alpha = 0.6,
                 position = "identity") +
  # Add mean lines for each time_group
  geom_vline(data = summary_stats_all, aes(xintercept = mean_length, color = time_group),
             linetype = "solid", linewidth = 1, show.legend = TRUE) + # Ensure legend for color is shown
  labs(x = NULL, y = "Count") +
  scale_fill_manual(name = "Survival Time", values = common_colors_length_fill) + # Distinct name for fill legend
  scale_color_manual(name = "Mean Length", # Distinct name for color legend (for mean lines)
                     values = common_colors_length_mean) +
  scale_x_continuous(breaks = x_breaks_val_length, limits = c(x_min_val_length, x_max_val_length)) +
  annotate("text", x = x_max_val_length, y = y_anno_all,
           label = "All samples", hjust = 1, vjust = 1,
           size = 5, fontface = "bold") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_blank(),
        legend.position = "none", 
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        panel.grid.major.x = element_line(color = "gray80", linetype = "dotted"),
        panel.grid.major.y = element_line(color = "gray80", linetype = "dotted"),
        plot.margin = unit(c(0, 5.5, 5.5, 5.5), "pt"))

# Sequenced samples only plot
p_sequenced_counts <- ggplot(sequenced_fish_modified, aes(x = total_length, fill = time_group)) +
  geom_histogram(binwidth = 2,
                 color = "black",
                 alpha = 0.6,
                 position = "identity") +
  # Add mean lines for each time_group
  geom_vline(data = summary_stats_sequenced, aes(xintercept = mean_length, color = time_group),
             linetype = "solid", linewidth = 1, show.legend = TRUE) + # Ensure legend for color is shown
  labs(x = "Total Length (mm)", y = "Count") +
  scale_fill_manual(name = "Survival Time", values = common_colors_length_fill) + # Distinct name
  scale_color_manual(name = "Mean Length", # Distinct name
                     values = common_colors_length_mean) +
  scale_x_continuous(breaks = x_breaks_val_length, limits = c(x_min_val_length, x_max_val_length)) +
  annotate("text", x = x_max_val_length, y = y_anno_sequenced,
           label = "sequenced", hjust = 1, vjust = 1,
           size = 5, fontface = "bold") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_blank(),
        legend.position = "none", 
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        panel.grid.major.x = element_line(color = "gray80", linetype = "dotted"),
        panel.grid.major.y = element_line(color = "gray80", linetype = "dotted"),
        plot.margin = unit(c(0, 5.5, 5.5, 5.5), "pt"))

# Combine Plots
combined_length_plots <- p_all_counts / p_sequenced_counts +
  plot_layout(guides = 'collect') & theme(legend.position = 'right') # Collect all legends and place them on the right

print(combined_length_plots)

# Save the combined plot
ggsave(paste0(save_dir, "Pmar_tfm_length_distributions_combined_filt_", plot_version, ".png"), combined_length_plots, dpi = 1000, height = 6, width = 8)
ggsave(paste0(save_dir, "Pmar_tfm_length_distributions_combined_filt_", plot_version, ".svg"), combined_length_plots, dpi = 1000, height = 6, width = 8)


#############################
# --- Plot 4: Per-Tank Length Distributions (SI) ---
#############################

# Get unique experiment names
experiment_names <- unique(all_fish_modified$experiment)

# Loop through each experiment
for (exp_name in experiment_names) {
  
  # Get unique tank names within the current experiment
  tanks_in_experiment <- all_fish_modified %>%
    filter(experiment == exp_name) %>%
    distinct(tank) %>%
    pull(tank)
  
  # Loop through each tank within the current experiment
  for (tank_name in tanks_in_experiment) {
    
    # Filter data for the current experiment and tank
    current_tank_data <- all_fish_modified %>%
      filter(experiment == exp_name, tank == tank_name)
    
    # Filter data for "sequenced" samples within the current experiment and tank
    sequenced_tank_data <- current_tank_data %>%
      filter(status == "sequenced")
    
    # Skip if no sequenced samples in this tank
    if (nrow(sequenced_tank_data) == 0) {
      message(paste("Skipping", exp_name, "tank", tank_name, "as no 'sequenced' samples were found."))
      next # Skip to the next iteration of the inner loop
    }
    
    # Calculate max counts for annotations for the current tank's data
    hist_data_all_tank <- hist(current_tank_data$total_length,
                               breaks = seq(x_min_val_length, x_max_val_length, by = 2), plot = FALSE)
    y_anno_all_tank <- max(hist_data_all_tank$counts) * 0.95
    
    hist_data_sequenced_tank <- hist(sequenced_tank_data$total_length,
                                     breaks = seq(x_min_val_length, x_max_val_length, by = 2), plot = FALSE)
    y_anno_sequenced_tank <- max(hist_data_sequenced_tank$counts) * 0.95
    
    # Calculate means for each time_group within the current tank's data
    summary_stats_all_tank <- current_tank_data %>%
      group_by(time_group) %>%
      summarise(
        mean_length = mean(total_length, na.rm = TRUE)
      )
    
    summary_stats_sequenced_tank <- sequenced_tank_data %>%
      group_by(time_group) %>%
      summarise(
        mean_length = mean(total_length, na.rm = TRUE)
      )
    
    # All Samples for current tank
    p_all_counts_tank <- ggplot(current_tank_data, aes(x = total_length, fill = time_group)) +
      geom_histogram(binwidth = 2,
                     color = "black",
                     alpha = 0.6,
                     position = "identity") +
      # Importantly, set show.legend = FALSE here for the mean line to avoid duplicate legend parts
      geom_vline(data = summary_stats_all_tank, aes(xintercept = mean_length, color = time_group),
                 linetype = "solid", linewidth = 1, show.legend = FALSE) + 
      labs(x = NULL, y = "Count") +
      # Also set show.legend = FALSE for the fill manual scale
      scale_fill_manual(name = "Survival Time", values = common_colors_length_fill, guide = "none") + 
      scale_color_manual(name = "Mean", values = common_colors_length_mean, guide = "none") + # Ensure no guide here either
      scale_x_continuous(breaks = x_breaks_val_length, limits = c(x_min_val_length, x_max_val_length)) +
      annotate("text", x = x_max_val_length, y = y_anno_all_tank,
               label = paste(exp_name, " ", tank_name, "\nall samples"), hjust = 1, vjust = 1,
               size = 4, fontface = "bold") +
      theme_bw() +
      theme(axis.text = element_text(size = 12),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.y = element_text(size = 14),
            plot.title = element_blank(),
            legend.position = "none", # Keep this as 'none'
            panel.grid.major.x = element_line(color = "gray80", linetype = "dotted"),
            panel.grid.major.y = element_line(color = "gray80", linetype = "dotted"),
            plot.margin = unit(c(5.5, 5.5, 0, 5.5), "pt"))
    
    # Sequenced Samples only for current tank
    p_sequenced_counts_tank <- ggplot(sequenced_tank_data, aes(x = total_length, fill = time_group)) +
      geom_histogram(binwidth = 2,
                     color = "black",
                     alpha = 0.6,
                     position = "identity") +
      # Here, explicitly allow legend for the color aesthetic (mean line)
      geom_vline(data = summary_stats_sequenced_tank, aes(xintercept = mean_length, color = time_group),
                 linetype = "solid", linewidth = 1, show.legend = TRUE) +
      labs(x = "Total Length (mm)", y = "Count") +
      # Explicitly allow legend for the fill aesthetic (histogram bars)
      scale_fill_manual(name = "Survival Time", values = common_colors_length_fill, guide = "legend") + 
      scale_color_manual(name = "Mean", values = common_colors_length_mean, guide = "legend") + # Explicitly allow guide for color
      scale_x_continuous(breaks = x_breaks_val_length, limits = c(x_min_val_length, x_max_val_length)) +
      annotate("text", x = x_max_val_length, y = y_anno_sequenced_tank,
               label = paste(exp_name, " ", tank_name, "\nsequenced"), hjust = 1, vjust = 1,
               size = 4, fontface = "bold") +
      theme_bw() +
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 14),
            plot.title = element_blank(),
            legend.position = "none", # Keep this as 'none' for individual plot
            legend.title = element_text(size = 12),
            legend.text = element_text(size = 10),
            panel.grid.major.x = element_line(color = "gray80", linetype = "dotted"),
            panel.grid.major.y = element_line(color = "gray80", linetype = "dotted"),
            plot.margin = unit(c(0, 5.5, 5.5, 5.5), "pt"))
    
    # Combine Plots using Patchwork for current tank
    # The key is to apply theme(legend.position = 'right') *after* plot_layout(guides = 'collect')
    combined_length_plots_tank <- p_all_counts_tank / p_sequenced_counts_tank +
      plot_layout(guides = 'collect') & theme(legend.position = 'right')
    
    print(combined_length_plots_tank)
    
    # Save the combined plot for the current experiment and tank
    ggsave(paste0(save_dir, "Pmar_tfm_length_distributions_", exp_name, "_", tank_name, "_with_stats_", plot_version, ".png"),
           combined_length_plots_tank, dpi = 1000, height = 6, width = 8)
    ggsave(paste0(save_dir, "Pmar_tfm_length_distributions_", exp_name, "_", tank_name, "_with_stats_", plot_version, ".svg"),
           combined_length_plots_tank, dpi = 1000, height = 6, width = 8)
  }
}
