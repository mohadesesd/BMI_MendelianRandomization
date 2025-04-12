# -----------------------------------------------------------
# Load Required Libraries
# -----------------------------------------------------------
library(data.table)  # For fast data manipulation
library(biomaRt)     # For querying Ensembl SNP data

# -----------------------------------------------------------
# 1. Read the sQTL and eQTL TSV Files
# -----------------------------------------------------------
# Update the file paths accordingly
sqtls <- fread("Filtered_sQTLs_KNOP1.tsv")
eqtls <- fread("Filtered_eQTLs_KNOP1.tsv")

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
# variant_dt now contains: variant, chr (as character), and pos.

# -----------------------------------------------------------
# 4. Query Ensembl (via biomaRt) to Retrieve rsIDs Using a ±100 bp Window
# -----------------------------------------------------------
# Connect to the Ensembl SNP mart.
mart <- useMart("ENSEMBL_MART_SNP",
                dataset = "hsapiens_snp",
                host = "https://www.ensembl.org")

# Set the query window (in base pairs)
window <- 100

# Prepare an empty list to store results for each chromosome
snp_list <- list()

# Loop over each unique chromosome in variant_dt:
for (ch in unique(variant_dt$chr)) {
  dt_ch <- variant_dt[chr == ch]
  pos_vec <- dt_ch$pos
  
  # Query Ensembl using filters "chr_name", "start", and "end"
  # and attributes: "refsnp_id", "chr_name", and "chrom_start".
  res <- tryCatch({
    getBM(
      attributes = c("refsnp_id", "chr_name", "chrom_start"),
      filters    = c("chr_name", "start", "end"),
      values     = list(
        chr_name = ch,
        start    = pos_vec - window,
        end      = pos_vec + window
      ),
      mart       = mart
    )
  }, error = function(e) {
    message("Error querying for chromosome ", ch, ": ", e)
    NULL
  })
  
  if (!is.null(res) && nrow(res) > 0) {
    snp_list[[ch]] <- as.data.table(res)
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

# Before merging, ensure that the join columns are the same type.
# Convert snp_mapping$chr_name to character (it might be integer).
snp_mapping[, chr_name := as.character(chr_name)]

# Rename the attribute "chrom_start" to "pos", if present.
setnames(snp_mapping, "chrom_start", "pos", skip_absent = TRUE)

# -----------------------------------------------------------
# 5. Merge the SNP Mapping with the Unique Variant Table
# -----------------------------------------------------------
# Merge variant_dt (with columns "chr", "pos", "variant") with snp_mapping (with "chr_name" and "pos")
variant_dt <- merge(variant_dt, snp_mapping, by.x = c("chr", "pos"),
                    by.y = c("chr_name", "pos"), all.x = TRUE)

# Rename the column "refsnp_id" to "SNP".
setnames(variant_dt, "refsnp_id", "SNP", skip_absent = TRUE)
# variant_dt now has columns: chr, pos, variant, and SNP.

# -----------------------------------------------------------
# 6. Merge the rsID Information Back into the sQTL and eQTL Data
# -----------------------------------------------------------
sqtls <- merge(sqtls, variant_dt[, .(variant, SNP)], by = "variant", all.x = TRUE)
eqtls <- merge(eqtls, variant_dt[, .(variant, SNP)], by = "variant", all.x = TRUE)

# -----------------------------------------------------------
# 7.  Save the Updated Tables as TSV Files
# -----------------------------------------------------------
fwrite(sqtls, file = "sQTL_with_SNP.tsv", sep = "\t")
fwrite(eqtls, file = "eQTL_with_SNP.tsv", sep = "\t")

cat("Finished processing. The updated sQTL and eQTL files have been saved as 'sQTL_with_SNP.tsv' and 'eQTL_with_SNP.tsv'.\n")
