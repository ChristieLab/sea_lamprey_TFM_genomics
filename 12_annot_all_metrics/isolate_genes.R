#isolate genes in outlier windows 
#NJCB

library(ggplot2)
library(qqman)
library(cowplot)
library(rtracklayer)
library(GenomicRanges)
library(writexl)
library(dplyr)

#install_github("stephenturner/qqman")


wd <- c("~/Purdue/Pmart_wgs/results/12_annot_all_metrics/")
setwd(wd)
list.files()


########################################################
#########        P. marinus early vs late      #########
#########           6.23 million SNPs           ########
########################################################

outlier_wins <- read.table("./outlier_windows.txt", sep = "\t", header = TRUE)

#Ok, now we will use the sea lamprey annotation to pull all genes that have been annotated that are within or overlap with these windows.
#Read in the GFF file
gff <- "~/Purdue/Pmart_wgs/reference/genomic.gff" # Replace with the actual path to your GFF file
lamprey_genes <- rtracklayer::import(gff)
lamprey_genes_filtered <- lamprey_genes[!grepl("^NW_", seqnames(lamprey_genes))]

#Make the windows a Granges object
windows_gr <- GRanges(
  seqnames = outlier_wins$chromosome,
  ranges = IRanges(start = outlier_wins$start, end = outlier_wins$end)
)

#check
seqlevels(windows_gr) <- as.character(seqlevels(windows_gr))
seqlevels(lamprey_genes) <- as.character(seqlevels(lamprey_genes))

common_seqlevels <- intersect(seqlevels(windows_gr), seqlevels(lamprey_genes))
windows_gr <- keepSeqlevels(windows_gr, common_seqlevels, pruning.mode = "coarse")
lamprey_genes <- keepSeqlevels(lamprey_genes, common_seqlevels, pruning.mode = "coarse")


#Find overlaps between the windows and the genes
overlaps <- GenomicRanges::findOverlaps(windows_gr, lamprey_genes)

#Extract the genes
overlapping_genes_info <- lamprey_genes[subjectHits(overlaps)]

# Convert the GRanges object to a data frame for easier filtering
overlapping_genes_df <- as.data.frame(overlapping_genes_info)

# Filter for rows where the 'type' column is equal to "gene"
gene_annotations <- overlapping_genes_df %>%
  filter(type == "gene")
mRNA_annotations <- overlapping_genes_df %>%
  filter(type == "mRNA")

#save for now
#write_xlsx(gene_annotations, path = "~/Purdue/Pmart_wgs/results/12_annot_all_metrics/outlier_windows_genes.xlsx")
#write_xlsx(mRNA_annotations, path = "~/Purdue/Pmart_wgs/results/12_annot_all_metrics/outlier_windows_mRNA.xlsx")

#clean up some single exon genes
gene_annotations_filt <- gene_annotations %>%
  filter(gene_annotations$width >= 200)
 # 
mRNA_annotations_filt <- mRNA_annotations %>%
  filter(mRNA_annotations$width >= 200)
 # 

#write_xlsx(gene_annotations_filt, path = "~/Purdue/Pmart_wgs/results/12_annot_all_metrics/outlier_windows_genes_min200bp.xlsx")
#write_xlsx(mRNA_annotations_filt, path = "~/Purdue/Pmart_wgs/results/12_annot_all_metrics/outlier_windows_mRNA_min200bp.xlsx")


#make list of parent RNAs
parents <- mRNA_annotations_filt$ID
#write.table(parents, file = "~/Purdue/Pmart_wgs/results/12_annot_all_metrics/parents.txt", sep = "\n", quote = FALSE, row.names = FALSE)


####move to blastp of the parents

#Get the top hits for each gene
######
final_data <- gene_annotations_filt %>%
  dplyr::select(ID)

# Read the BLAST results and select the best hit for each mRNA
top_hits <- read.csv("~/Purdue/Pmart_wgs/results/12_annot_all_metrics/outlier_window_mRNA_blastp_results_concat/best_hits/Pmart_outlier_window_blast_results_combined_best_hits.csv", header = TRUE)

best_mrna_hits <- top_hits %>%
  arrange(evalue, desc(bitscore)) %>%
  group_by(qseqid) %>%
  slice(1) %>%
  ungroup()

# Join the best mRNA hits with the mRNA annotations to get location and parent info
best_mrna_hits_with_genes <- best_mrna_hits %>%
  dplyr::left_join(mRNA_annotations_filt, by = c("qseqid" = "ID"))

# Convert the Parent column from a list to a character vector to allow for the join
best_mrna_hits_with_genes$Parent <- sapply(best_mrna_hits_with_genes$Parent, "[[", 1)

# Now, join the best hits back to the original filtered gene list using the parent gene ID
best_gene_hits_df <- final_data %>%
  dplyr::left_join(best_mrna_hits_with_genes, by = c("ID" = "Parent"))

# If a gene had multiple mRNAs with hits, this will select the best one based on e-value/bitscore
best_gene_hits_df <- best_gene_hits_df %>%
  group_by(ID) %>%
  arrange(evalue, desc(bitscore)) %>%
  slice(1) %>%
  ungroup()

# Define the new column order 
final_column_order <- c(
  "ID", 
  "qseqid",
  "seqnames",
  "start",
  "end",
  "width",
  "strand",
  "stitle",
  "sseqid",
  "query_coverage",
  "evalue",
  "pident",
  "bitscore",
  "length",
  "mismatch",
  "gapopen",
  "qstart",
  "qend",
  "sstart",
  "send",
  "aligned_query_length",
  "query_length"
)

# Reorder the columns 
final_blast_results_with_parent_df <- best_gene_hits_df %>%
  dplyr::select(dplyr::all_of(final_column_order)) %>%
  dplyr::rename(Parent = ID) %>%
  dplyr::mutate(Parent = as.character(Parent)) %>%
  dplyr::arrange(seqnames, start)

# Now, write the reordered and cleaned data frame to an Excel file
write_xlsx(final_blast_results_with_parent_df, path = "~/Purdue/Pmart_wgs/results/12_annot_all_metrics/outlier_window_FINAL_best_gene_hit_blast_results.xlsx")



