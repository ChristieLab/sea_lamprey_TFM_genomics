library(ggplot2)
library(dplyr)
library(ggrepel)

# Read in data
setwd("~/Purdue/Pmart_wgs/results/14_karyotpye/")
outlier_chr <- read.csv("./Pmar_chrom_size_outlier.csv", header = TRUE)
outlier_win <- read.table("./outlier_windows.txt", header = TRUE)

## 1. index y
outlier_chr$y_index <- rev(seq_len(nrow(outlier_chr)))

## Map outlier windows 
outlier_win$y_val <- outlier_chr$y_index[match(outlier_win$chromosome, outlier_chr$chr)]

###############################
## Create the Karyotype Plot ##
###############################
p <- ggplot() +
  # Draw chromosome bars
  geom_rect(data = outlier_chr,
            aes(xmin = 0, xmax = length, 
                ymin = y_index - 0.4, ymax = y_index + 0.4),
            fill = "#0D98BA", color = "grey40", linewidth = 0.2, alpha = 0.2) +
  # Draw outlier highlight regions
  geom_rect(data = outlier_win,
            aes(xmin = start, xmax = end, 
                ymin = y_val - 0.35, ymax = y_val + 0.35),
            fill = "#01386A", alpha = 1.0) +
  # Set Y-axis labels to use the 'num' column
  scale_y_continuous(breaks = outlier_chr$y_index, 
                     labels = outlier_chr$num) +
  # Format X-axis for Megabases
  scale_x_continuous(labels = function(x) paste0(x / 1e6, " Mb")) +
  labs(x = "Genomic Position", 
       y = "Outlier Chromosome") +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.background = element_rect(fill = "white", color = NA))

p

ggsave(paste("~/Purdue/Pmart_wgs/figures_and_tables/14_karyotpye/Fig_5_karyotype_panel_v1.png"), p, dpi = 1000, height = 4, width = 6.5)

#################################################################
###lets add in some ways to label certain categories of genes. ##
#################################################################

out_cat <- read.csv("./outlier_genes_category_II_short.csv", header = TRUE)
top4 <- read.csv("./Top4outliergene_Category.csv", header = TRUE)

# Process out_cat coordinates
out_cat$midpoint <- (out_cat$start + out_cat$end) / 2
out_cat$y_val <- outlier_chr$y_index[match(out_cat$chrom, outlier_chr$chr)]
out_cat$category <- as.factor(out_cat$category)

# Process top4 coordinates for lollipops and labels
top4$midpoint <- (top4$start + top4$end) / 2
top4$y_val <- outlier_chr$y_index[match(top4$chromosome, outlier_chr$chr)]
top4$Category <- as.factor(top4$Category)

cat_labels <- c("1" = "Mitochondrial function & energy metabolism",
                "2" = "Detoxification & transport",
                "3" = "Cellular stress & apoptosis") 

bright_palette <- c(
  "1" = "#FFD700", 
  "2" = "#FF00FF", 
  "3" = "#FF4500" 
)

p3 <- ggplot() +
  # Draw chromosome bars
  geom_rect(data = outlier_chr,
            aes(xmin = 0, xmax = length, 
                ymin = y_index - 0.4, ymax = y_index + 0.4),
            fill = "#0D98BA", color = "grey40", linewidth = 0.3, alpha = 0.2) +
  
  # Draw outlier highlight regions
  geom_rect(data = outlier_win,
            aes(xmin = start, xmax = end, 
                ymin = y_val - 0.35, ymax = y_val + 0.35),
            fill = "#01386A", alpha = 0.5) +
  
  # Add gene category segments (Barcode Style)
  geom_segment(data = out_cat, 
               aes(x = midpoint, xend = midpoint, 
                   y = y_val - 0.35, yend = y_val + 0.35, 
                   color = category),
               linewidth = 2) + 
  
  # Add Lollipops for top4 genes
  geom_segment(data = top4,
               aes(x = midpoint, xend = midpoint, 
                   y = y_val + 0.35, yend = y_val + 0.6),
               color = "grey30", linewidth = 0.2) +
  geom_point(data = top4,
             aes(x = midpoint, y = y_val + 0.55, fill = Category),
             shape = 21, size = 3, color = "black", stroke = 0.5) +
  
  # Add labels for top4 genes
  geom_text_repel(data = top4,
                  aes(x = midpoint, y = y_val + 0.6, label = Protein_abbreviation),
                  size = 3, fontface = "italic", nudge_y = 0.3) +
  
  # Aesthetics and Labels
  scale_color_manual(values = bright_palette, labels = cat_labels) +
  scale_fill_manual(values = bright_palette, labels = cat_labels) +
  scale_y_continuous(breaks = outlier_chr$y_index, 
                     labels = outlier_chr$num) +
  scale_x_continuous(labels = function(x) paste0(x / 1e6, " Mb")) +
  
  labs(x = "Genomic Position", 
       y = "Outlier Chromosome",
       color = "Functional Category",
       fill = "Functional Category") +
  
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        legend.position = c(0.80, 0.25))

p3

ggsave(paste("~/Purdue/Pmart_wgs/figures_and_tables/14_karyotpye/Fig_5_karyotype_panel_v1.4.png"), p3, dpi = 1000, height = 6, width = 8)
ggsave(paste("~/Purdue/Pmart_wgs/figures_and_tables/14_karyotpye/Fig_5_karyotype_panel_v1.4.svg"), p3, dpi = 1000, height = 6, width = 8)
