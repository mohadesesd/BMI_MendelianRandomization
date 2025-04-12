#!/usr/bin/env Rscript
library(dplyr)
library(readr)
library(stringr)

#-----------------------------------------------------
# Hard-coded file locations:
#-----------------------------------------------------
mr_dir      <- "LCL_MR/Raw_MR/"              # Directory containing MR_Chr*.rds files
mapping_rds <- "mapping_table.rds"                # Pre-saved mapping table RDS file
output_dir  <- "LCL_MR/Cleaned_Raw_MR/"         # Directory to save cleaned RDS files

# Create the output directory if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

#-----------------------------------------------------
# Step 1: Load and combine MR results from all RDS files
#-----------------------------------------------------
mr_files <- list.files(path = mr_dir, pattern = "^Raw_MR_Chr.*\\.rds$", full.names = TRUE)
if(length(mr_files) == 0) {
  stop("No MR_Chr*.rds files found in the specified directory.")
}

# Read and combine all MR files
mr_list <- lapply(mr_files, readRDS)
mr_data <- bind_rows(mr_list)

# Check for required columns
if(!("id.exposure" %in% colnames(mr_data))) {
  stop("Combined MR results do not contain the 'id.exposure' column.")
}
if(!("exposure" %in% colnames(mr_data))) {
  stop("Combined MR results do not contain the 'exposure' column.")
}

#-----------------------------------------------------
# Step 2: Extract chromosome and gene_id from the sQTL identifier
#-----------------------------------------------------
# Assumes the identifier format: "chrX:start:end:clu_xxxx:gene_id.version"
mr_data <- mr_data %>%
  mutate(
    chr_from_id = str_split(id.exposure, ":", simplify = TRUE)[,1],
    raw_gene_id = str_split(id.exposure, ":", simplify = TRUE)[,5],
    gene_id     = sub("\\..*", "", raw_gene_id)  # Remove version (e.g., ".14")
  )

#-----------------------------------------------------
# Step 3: Load the pre-saved mapping table
#-----------------------------------------------------
mapping_table <- readRDS(mapping_rds)
# The mapping_table is assumed to have columns: chr, gene_id, gene_name, and biotype

#-----------------------------------------------------
# Step 4: Join the mapping table to add gene mapping info to MR results
#-----------------------------------------------------
mr_data <- mr_data %>%
  left_join(mapping_table, by = "gene_id")

#-----------------------------------------------------
# Step 5: Determine the final chromosome for each instrument
#-----------------------------------------------------
# Use the chromosome from the mapping table (column "chr") if available; otherwise, use chr_from_id.
mr_data <- mr_data %>%
  mutate(chr_final = if_else(!is.na(chr), chr, chr_from_id))

#-----------------------------------------------------
# Step 6: For exposures with multiple instruments, filter out those not on the dominant chromosome
#-----------------------------------------------------
# For each exposure, compute the most common (dominant) chr_final value and retain only instruments on that chromosome.
mr_data <- mr_data %>%
  group_by(exposure) %>%
  mutate(dominant_chr = names(which.max(table(chr_final)))) %>%
  ungroup() %>%
  filter(chr_final == dominant_chr)

#-----------------------------------------------------
# Step 7: Split the cleaned data by chromosome and save each subset as an RDS file
#-----------------------------------------------------
unique_chrs <- unique(mr_data$chr_final)
for(ch in unique_chrs) {
  subset_data <- mr_data %>% filter(chr_final == ch)
  output_file <- file.path(output_dir, paste0("cleaned_MR_", ch, ".rds"))
  saveRDS(subset_data, file = output_file)
  cat("Saved cleaned data for chromosome", ch, "to", output_file, "\n")
}
