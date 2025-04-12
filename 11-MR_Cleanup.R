#!/usr/bin/env Rscript
library(dplyr)
library(readr)

#---------------------------
# File paths (modify as needed)
#---------------------------
gtf_file   <- "Homo_sapiens.GRCh38.112.gene.gtf"  # Your GTF file
output_rds <- "mapping_table.rds"                # Output RDS file for the mapping table

#---------------------------
# Step 1: Read the GTF file
#---------------------------
# Read the file as tab-delimited; it has no header and comment lines start with "#"
gtf <- read.delim(gtf_file, header = FALSE, comment.char = "#", stringsAsFactors = FALSE)
colnames(gtf) <- c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

#---------------------------
# Step 2: Filter for gene features
#---------------------------
gtf_genes <- gtf %>% filter(feature == "gene")

#---------------------------
# Step 3: Define a function to parse the attribute column
#---------------------------
parse_attributes <- function(attr) {
  # Split by semicolon, remove any empty strings and trim whitespace
  pairs <- unlist(strsplit(attr, ";"))
  pairs <- trimws(pairs)
  pairs <- pairs[pairs != ""]
  
  # Initialize an empty list to hold key-value pairs
  res <- list()
  for(p in pairs) {
    # Split by whitespace; assume format is: key "value"
    kv <- unlist(strsplit(p, " ", fixed = TRUE))
    if(length(kv) >= 2) {
      key <- kv[1]
      value <- kv[2]
      # Remove any quotes from the value
      value <- gsub('"', '', value)
      res[[key]] <- value
    }
  }
  return(res)
}

#---------------------------
# Step 4: Apply the parser to each row's attribute
#---------------------------
parsed_attrs <- lapply(gtf_genes$attribute, parse_attributes)

# Extract gene_id, gene_name, and gene_biotype into vectors
gene_id_vec   <- sapply(parsed_attrs, function(x) if("gene_id" %in% names(x)) x[["gene_id"]] else NA)
gene_name_vec <- sapply(parsed_attrs, function(x) if("gene_name" %in% names(x)) x[["gene_name"]] else NA)
biotype_vec   <- sapply(parsed_attrs, function(x) if("gene_biotype" %in% names(x)) x[["gene_biotype"]] else NA)

#---------------------------
# Step 5: Create the mapping table data frame
#---------------------------
mapping_table <- data.frame(
  chr = gtf_genes$chr,
  gene_id = gene_id_vec,
  gene_name = gene_name_vec,
  biotype = biotype_vec,
  stringsAsFactors = FALSE
)

# Standardize the chromosome names: if not starting with "chr", add the prefix.
mapping_table <- mapping_table %>%
  mutate(chr = if_else(grepl("^chr", chr), chr, paste0("chr", chr)))

#---------------------------
# Step 6: Save the mapping table as an RDS file
#---------------------------
saveRDS(mapping_table, file = output_rds)
cat("Mapping table saved to", output_rds, "\n")
