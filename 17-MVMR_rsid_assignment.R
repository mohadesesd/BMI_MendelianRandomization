# -----------------------------------------------------------
# Load Required Libraries
# -----------------------------------------------------------
library(curl)         # For making HTTP requests
library(data.table)   # For fast data manipulation
library(jsonlite)     # For parsing JSON responses

# -----------------------------------------------------------
# 1. Read the sQTL and eQTL TSV Files
# -----------------------------------------------------------
# Update the file paths accordingly
sqtls <- fread("Filtered_sQTLs_FKBP11.tsv")
eqtls <- fread("Filtered_eQTLs_FKBP11.tsv")

# -----------------------------------------------------------
# 2. Create a "variant" Column as "chr:pos"
# -----------------------------------------------------------
# Assumes each file already has columns "chr" and "pos"
sqtls[, variant := paste0(chr, ":", pos)]
eqtls[, variant := paste0(chr, ":", pos)]

# -----------------------------------------------------------
# 3. Get the Unique Variants Across Both Datasets
# -----------------------------------------------------------
all_variants <- unique(c(sqtls$variant, eqtls$variant))
variant_dt <- data.table(variant = all_variants)
# Split the 'variant' string into "chr" and "pos" columns.
variant_dt[, c("chr", "pos") := tstrsplit(variant, ":", fixed = TRUE)]
variant_dt[, pos := as.numeric(pos)]
# For the Ensembl query, remove the "chr" prefix.
variant_dt[, chr := sub("^chr", "", chr)]

# -----------------------------------------------------------
# 4. Query Ensembl API (Overlap Endpoint) Using curl to Retrieve SNP IDs (rsIDs)
# -----------------------------------------------------------
query_ensembl_curl <- function(chromosome, pos_vec, window = 100, retries = 3, wait_time = 5) {
  # Use the overlap endpoint which returns all features overlapping the region.
  ensembl_url <- "https://rest.ensembl.org/overlap/region/human/"
  
  # Prepare the query string.
  query <- paste0(chromosome, ":", min(pos_vec) - window, "-", max(pos_vec) + window, "?feature=variation")
  
  # Debug: Print the URL to check if it's correct.
  full_url <- paste0(ensembl_url, query)
  message("Querying URL: ", full_url)
  
  # Prepare the curl handle with headers.
  h <- new_handle(timeout = 300)
  handle_setheaders(h,
                    "Content-Type" = "application/json",
                    "Accept" = "application/json")
  
  # Retry mechanism.
  attempt <- 1
  while (attempt <= retries) {
    res <- tryCatch({
      # Fetch the URL using curl.
      r <- curl_fetch_memory(full_url, handle = h)
      raw_response <- rawToChar(r$content)
      
      # Check if response is HTML (error page)
      if (grepl("<html>", raw_response, ignore.case = TRUE)) {
        stop("Received an HTML error page. Possibly incorrect request.")
      }
      
      # Parse JSON response.
      snp_data <- fromJSON(raw_response)
      return(snp_data)
    }, error = function(e) {
      message("Error querying for chromosome ", chromosome, ": ", e)
      return(NULL)
    })
    
    if (!is.null(res) && length(res) > 0) {
      return(res)
    }
    
    message(paste("Query failed. Retrying... (Attempt", attempt, "of", retries, ")"))
    Sys.sleep(wait_time)
    attempt <- attempt + 1
  }
  return(NULL)
}

# Prepare an empty list to store results for each chromosome.
snp_list <- list()

# Loop over each unique chromosome in variant_dt:
for (ch in unique(variant_dt$chr)) {
  dt_ch <- variant_dt[chr == ch]
  pos_vec <- dt_ch$pos
  
  # Query Ensembl with retries.
  res <- query_ensembl_curl(ch, pos_vec)
  
  if (!is.null(res) && length(res) > 0) {
    # 'res' should be a JSON array of objects.
    # We expect each object to have: "id" (rsID), "seq_region_name", and "start".
    # Create a data.table with:
    #   - variant: constructed as "seq_region_name_start"
    #   - SNP: the id field (which is typically something like "rsXXXXX")
    snp_dt <- data.table(
      variant = paste0(res$seq_region_name, "_", res$start),
      SNP = res$id
    )
    snp_list[[ch]] <- snp_dt
  } else {
    message("No mapping returned for chromosome ", ch)
  }
}

# Combine the results for all chromosomes.
if (length(snp_list) == 0) {
  message("No SNP mapping results obtained from Ensembl. Creating empty mapping.")
  snp_mapping <- data.table(chr_name = character(),
                            pos = numeric(),
                            refsnp_id = character())
} else {
  snp_mapping <- rbindlist(snp_list, fill = TRUE)
}

# For the merge, we need to have join columns.
# In our constructed snp_dt, the column 'variant' is created as "seq_region_name_start".
# For clarity, we add two columns: 'chr_name' and 'pos' extracted from that.
if (nrow(snp_mapping) > 0) {
  snp_mapping[, chr_name := as.character(sapply(variant, function(x) strsplit(x, "_")[[1]][1]))]
  snp_mapping[, pos := as.numeric(sapply(variant, function(x) strsplit(x, "_")[[1]][2]))]
}

message("Columns in snp_mapping: ", paste(names(snp_mapping), collapse = ", "))

# -----------------------------------------------------------
# 5. Merge the SNP Mapping with the Unique Variant Table Using .EACHI for Group-wise Join
# -----------------------------------------------------------
variant_dt <- variant_dt[
  snp_mapping, 
  on = .(chr = chr_name, pos = pos), 
  .(variant, SNP = SNP), 
  by = .EACHI
]

# -----------------------------------------------------------
# 6. Merge the rsID Information Back into the sQTL and eQTL Data
# -----------------------------------------------------------
sqtls <- merge(sqtls, variant_dt[, .(variant, SNP)], by = "variant", all.x = TRUE)
eqtls <- merge(eqtls, variant_dt[, .(variant, SNP)], by = "variant", all.x = TRUE)

# -----------------------------------------------------------
# 7. Save the Updated Tables as TSV Files
# -----------------------------------------------------------
fwrite(sqtls, file = "sQTL_with_SNP.tsv", sep = "\t")
fwrite(eqtls, file = "eQTL_with_SNP.tsv", sep = "\t")

cat("Finished processing. The updated sQTL and eQTL files have been saved as 'sQTL_with_SNP.tsv' and 'eQTL_with_SNP.tsv'.\n")
