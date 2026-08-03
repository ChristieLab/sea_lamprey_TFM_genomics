#Filter blast outputs 

wd <- "/scratch/bell/nbackens/Pmart_wgs/vcf_investigation/annot_all_metrics/outlier_window_mRNA_blastp_results_concat/"
setwd(wd)
# Load necessary packages
library(dplyr)

# --- Configuration ---
blast_output_directory <- "/scratch/bell/nbackens/Pmart_wgs/vcf_investigation/annot_all_metrics/outlier_window_mRNA_blastp_results_concat/" 
fasta_sequences_directory <- "/scratch/bell/nbackens/Pmart_wgs/vcf_investigation/annot_all_metrics/outlier_window_proteins_split/" 
output_directory <- "/scratch/bell/nbackens/Pmart_wgs/vcf_investigation/annot_all_metrics/outlier_window_mRNA_blastp_results_concat/best_hits" 
combined_output_file <- file.path(output_directory, "Pmart_outlier_window_blast_results_combined_best_hits.csv")


# Create the output directory if it doesn't exist
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}

# Get a list of all .out files in the specified directory
blast_files <- list.files(blast_output_directory, pattern = "\\.out$", full.names = TRUE)

# Function to read FASTA and get lengths (without Biostrings)
get_fasta_lengths <- function(fasta_file) {
  sequences <- readLines(fasta_file)
  header_lines <- grep("^>", sequences)
  length_df <- data.frame(qseqid = character(), query_length = integer(),
                          stringsAsFactors = FALSE)
  
  for (i in seq_along(header_lines)) {
    start_line <- header_lines[i] + 1
    end_line <- ifelse(i < length(header_lines), header_lines[i + 1] - 1, length(sequences))
    sequence <- paste(sequences[start_line:end_line], collapse = "")
    header <- gsub(">", "", sequences[start_line - 1])
    # Extract the sequence ID (assuming it's the first word before a space)
    qseqid <- strsplit(header, " ")[[1]][1]
    length_df <- rbind(length_df, data.frame(qseqid = qseqid, query_length = nchar(sequence),
                                             stringsAsFactors = FALSE))
  }
  return(length_df)
}

# Function to process a single BLAST output file
process_blast_file <- function(blast_file) {
  # 1. Read the BLAST output file
  column_names <- c("qseqid", "sseqid", "pident", "length", "mismatch",
                    "gapopen", "qstart", "qend", "sstart", "send",
                    "evalue", "bitscore", "sseqid", "stitle")
  blast_data <- read.table(blast_file, header = FALSE, sep = "\t",
                           comment.char = "#", quote = "",
                           col.names = column_names,
                           stringsAsFactors = FALSE)
  
  if (nrow(blast_data) == 0) {
    cat("Warning: BLAST output file", basename(blast_file), "is empty. Skipping.\n")
    return(NULL)
  }
  
  # 2. Get the corresponding FASTA file path
  base_name <- gsub("\\.out$", "", basename(blast_file))
  fasta_file_name <- paste0(base_name, ".fa")
  fasta_file_path <- file.path(fasta_sequences_directory, fasta_file_name)
  
  # 3. Read the FASTA file and get query lengths (using our function)
  if (!file.exists(fasta_file_path)) {
    cat("Error: Corresponding FASTA file", fasta_file_path, "not found. Skipping", basename(blast_file), "\n")
    return(NULL)
  }
  
  query_lengths <- get_fasta_lengths(fasta_file_path)
  
  # 4. Calculate query coverage
  blast_data$aligned_query_length <- abs(blast_data$qend - blast_data$qstart) + 1
  blast_data <- merge(blast_data, query_lengths, by = "qseqid", all.x = TRUE)
  if (any(is.na(blast_data$query_length))) {
    warning("Some query IDs in", basename(blast_file), "not found in FASTA. Query coverage for these will be NA.")
  }
  blast_data$query_coverage <- (blast_data$aligned_query_length / blast_data$query_length) * 100
  
  # 5. Filter out unwanted stitles
  unwanted_strings <- c("uncharacterized", "hypothetical", "unnamed")
  blast_data_filtered <- blast_data %>%
    filter(!grepl(paste(unwanted_strings, collapse = "|"), stitle, ignore.case = TRUE))
  
  # 6. Select the top 5 hits
  top_5_hits <- blast_data_filtered %>%
    group_by(qseqid) %>%
    arrange(desc(query_coverage), desc(pident), evalue, .keep_all = TRUE) %>%
    slice_head(n = 5) %>%
    ungroup()
  
  # 7. Write the top 5 hits to an individual file
  output_file_name <- paste0(base_name, "_best_hits.tsv")
  output_file_path <- file.path(output_directory, output_file_name)
  write.table(top_5_hits, file = output_file_path, sep = "\t",
              row.names = FALSE, col.names = TRUE, quote = FALSE)
  cat("Top 5 hits for", basename(blast_file), "written to", output_file_path, "\n")
  
  return(top_5_hits)
}

# Process each BLAST file in the directory
all_best_hits <- lapply(blast_files, process_blast_file)

# Combine all the best hits into a single data frame
combined_best_hits_df <- do.call(rbind, all_best_hits)
combined_best_hits_df$stitle <- gsub(",", "", combined_best_hits_df$stitle)
# Write the combined best hits to a single file
write.csv(combined_best_hits_df, file = combined_output_file, quote = FALSE, col.names = TRUE)



