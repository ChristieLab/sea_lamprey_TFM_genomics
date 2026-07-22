#merge genes from annotation  
#NJCB

library(tidyverse)

wd <- c("~/Purdue/Pmart_wgs/results/12_annot_all_metrics/")
setwd(wd)
list.files()

annot <- readxl::read_excel("./outlier_window_FINAL_best_gene_hit_blast_results.xlsx")
symbol <- read.table("./gene_tsv_uniprot_id_mapping_v2/idmapping_2025_11_06.tsv", sep = "\t", quote = "", header = TRUE)
colnames(symbol)
colnames(symbol) <- c("sseqid", "Entry", "Reviewed", "Entry.Name", "Protein.names", "Gene.Names", "Organism", "Length")

merged_annot <- annot %>%
  left_join(symbol, by = "sseqid")

final_report <- merged_annot %>%
  select(seqnames, start, end, Gene.Names, stitle, sseqid)


colnames(final_report) <- c("chrom", "start", "end", "gene_symbol", "protein", "protein_accession")

writexl::write_xlsx(final_report, "./outlier_window_FINAL_best_gene_hit_blast_results_merged_with_uniprot_symbol.xlsx")
