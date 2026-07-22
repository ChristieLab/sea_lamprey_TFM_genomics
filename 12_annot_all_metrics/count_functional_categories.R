#sort categories of genes in outlier windows
#NJCB

library(readxl)
library(writexl)
library(tidyverse)

wd <- "~/Purdue/Pmart_wgs/results/12_annot_all_metrics/"
setwd(wd)

gene <- read_excel("./outlier_window_FINAL_best_gene_hit_blast_results_merged_with_uniprot_symbol_manual_categories_II.xlsx")

colnames(gene)

#List the top 50 terms found
top_functional_terms <- gene %>%
  select(functional_category, functional_category_2) %>%
  pivot_longer(cols = everything(), names_to = "source", values_to = "term") %>%
  filter(!is.na(term)) %>%  # Remove empty/NA cells
  count(term, sort = TRUE) %>%
  slice_max(n, n = 50)

top_functional_terms2 <- gene %>%
  select(functional_category, functional_category_2) %>%
  pivot_longer(cols = everything(), names_to = "source", values_to = "term") %>%
  filter(!is.na(term)) %>%
  # Standardize all terms to lowercase
  mutate(term = tolower(term)) %>% 
  count(term, sort = TRUE) %>%
  slice_max(n, n = 50)

# View the result
print(top_functional_terms2, n = 50)

write_xlsx(top_functional_terms2, "~/Purdue/Pmart_wgs/results/12_annot_all_metrics/outlier_window_FINAL_best_gene_hit_blast_results_merged_with_uniprot_symbol_manual_categories_II_TOP_CATEGORY_COUNTS.xlsx", col_names = TRUE)


mito_genes_broad <- gene %>%
  filter(
    str_detect(functional_category, regex("mitochon", ignore_case = TRUE)) |
      str_detect(functional_category_2, regex("mitochon", ignore_case = TRUE))
  )

# Count the unique gene symbols
total_unique_mito <- mito_genes_broad %>% 
  distinct(gene_symbol) %>% 
  nrow()

