# -----------------------------------------------------------
# This script pre-processes outcome data from the main outcome file and proxied outcome data,
# reads processed QTL files (eQTL_with_SNP.tsv and sQTL_with_SNP.tsv) for exposures,
# deduplicates data by SNP to avoid explosive cartesian joins,
# merges the exposures with outcome associations,
# reads and harmonizes the LD matrix,
# regularizes the LD matrix to avoid singularity,
# creates an MR input object using the MendelianRandomization package,
# and runs multivariable IVW and MR-Egger analyses.
# -----------------------------------------------------------

# Load required libraries
library(data.table)            # Fast data manipulation
library(dplyr)                 # Data wrangling
library(MendelianRandomization) # For MR analyses

# -------------------------------
# 1. Pre-process Outcome Data (Using variant_id as rsID)
# -------------------------------
i <- 16  # Chromosome of interest

# Read the main outcome associations file
df_outcome <- fread("GCST90002409_buildGRCh38.tsv")

# Filter for chromosome 16
GCST2409 <- df_outcome %>% filter(sprintf("chr%s", i) == chromosome)

# Convert allele columns to uppercase
GCST2409$effect_allele <- toupper(GCST2409$effect_allele)
GCST2409$other_allele  <- toupper(GCST2409$other_allele)
cat("Outcome data (main) after allele conversion (first few rows):\n")
print(head(GCST2409))

# Use the provided variant_id (which are rsIDs) as the SNP identifier
GCST2409[, SNP := variant_id]

# Read proxied outcome data from an RDS file
proxied <- readRDS('./Adipose_vis_MR/outcomeToMerge.rds')
proxied <- proxied %>% select(-SNP_Original)
proxied <- proxied %>% filter(sprintf("%s", i) == chromosome)
setDT(GCST2409)
setDT(proxied)

# Diagnose column names
cat("Columns in main outcome data:\n")
print(names(GCST2409))
cat("Columns in proxied outcome data:\n")
print(names(proxied))

# Determine common columns
common_cols <- intersect(names(GCST2409), names(proxied))
cat("Common outcome columns:\n")
print(common_cols)

# Bind using only common columns
GCST2409_combined <- rbindlist(list(GCST2409[, common_cols, with = FALSE],
                                    proxied[, common_cols, with = FALSE]), fill = TRUE)
cat("Combined outcome data dimensions:", dim(GCST2409_combined), "\n")

# Remove indels: keep rows where both alleles are A, T, C, or G
GCST2409_combined <- GCST2409_combined[
  other_allele %in% c("A", "T", "C", "G") & 
  effect_allele %in% c("A", "T", "C", "G"), 
]

# Select and reorder columns of interest
desired_cols <- c("variant_id", "SNP", "effect_allele", "other_allele", "p_value", "beta", "standard_error", "TotalSampleSize")
existing_desired <- intersect(desired_cols, names(GCST2409_combined))
GCST2409_Renamed <- GCST2409_combined[, ..existing_desired]

# Rename columns to standard names
new_names <- c("variant_id", "SNP", "effect_allele.outcome", "other_allele.outcome", 
               "pval.outcome", "beta.outcome", "se.outcome", "samplesize.outcome")
new_names <- new_names[1:length(existing_desired)]
setnames(GCST2409_Renamed, old = existing_desired, new = new_names)

# Remove palindromic variants
GCST2409_Renamed <- GCST2409_Renamed[
  !((effect_allele.outcome == "A" & other_allele.outcome == "T") |
    (effect_allele.outcome == "T" & other_allele.outcome == "A") |
    (effect_allele.outcome == "C" & other_allele.outcome == "G") |
    (effect_allele.outcome == "G" & other_allele.outcome == "C")), 
]

# Add outcome labels
GCST2409_Renamed$outcome <- "BMI"
GCST2409_Renamed$id.outcome <- "BMI"
cat("Pre-processed outcome data (first few rows):\n")
print(head(GCST2409_Renamed))

# Deduplicate outcome data by SNP
GCST2409_Renamed <- unique(GCST2409_Renamed, by = "SNP")
cat("Outcome data after deduplication (first few rows):\n")
print(head(GCST2409_Renamed))

# -------------------------------
# 2. Read Processed QTL Files (Exposures)
# -------------------------------
eqtl_data <- fread("Final_IVs_eQTL.tsv")
sqtl_data <- fread("Final_IVs_sQTL.tsv")
eqtl_data <- eqtl_data[!is.na(SNP) & SNP != ""]
sqtl_data <- sqtl_data[!is.na(SNP) & SNP != ""]

# Deduplicate by SNP
eqtl_data <- unique(eqtl_data, by = "SNP")
sqtl_data <- unique(sqtl_data, by = "SNP")

# -------------------------------
# 3. Merge QTL Data to Form Exposure Estimates
# -------------------------------
merged_exposures <- merge(eqtl_data, sqtl_data, by = "SNP", 
                          suffixes = c("_expr", "_splice"), all = TRUE)
merged_exposures <- unique(merged_exposures, by = "SNP")

# -------------------------------
# 4. Merge Exposure Data with Outcome Data
# -------------------------------
merged_data <- merge(merged_exposures, GCST2409_Renamed, by = "SNP", 
                     all.x = TRUE, allow.cartesian = TRUE)
cat("First few rows of merged exposure + outcome data:\n")
print(head(merged_data))
# Filter out rows with missing outcome associations
merged_data <- merged_data[!is.na(beta.outcome) & !is.na(se.outcome)]
cat("Merged data after filtering for outcome associations (first few rows):\n")
print(head(merged_data))

# -------------------------------
# 5. Read the LD Matrix and Harmonize SNP IDs
# -------------------------------
ld_matrix_result <- read.table("LD_matrix_result.tsv", header = TRUE, sep = "\t", 
                               row.names = 1, stringsAsFactors = FALSE, fill = TRUE, check.names = FALSE)
ld_matrix_result <- as.matrix(ld_matrix_result)
cat("Original LD matrix dimensions:", dim(ld_matrix_result), "\n")
if(nrow(ld_matrix_result) != ncol(ld_matrix_result)){
  cat("LD matrix is not square; subsetting columns to match rows.\n")
  ld_matrix_result <- ld_matrix_result[, 1:nrow(ld_matrix_result)]
  cat("New LD matrix dimensions:", dim(ld_matrix_result), "\n")
}
ld_snps <- rownames(ld_matrix_result)
cat("Original SNP IDs from LD matrix:\n")
print(ld_snps)

# Harmonize SNP IDs: force merged_data SNP IDs to include the "rs" prefix
merged_data[, SNP := ifelse(grepl("^rs", SNP), SNP, paste0("rs", SNP))]
all_IVs <- merged_data$SNP
all_IVs_adj <- ifelse(grepl("^rs", all_IVs), all_IVs, paste0("rs", all_IVs))
ld_snps_adj <- ifelse(grepl("^rs", ld_snps), ld_snps, paste0("rs", ld_snps))
cat("Adjusted SNP IDs from merged data:\n")
print(all_IVs_adj)
cat("Adjusted SNP IDs from LD matrix:\n")
print(ld_snps_adj)

rownames(ld_matrix_result) <- ld_snps_adj
colnames(ld_matrix_result) <- ld_snps_adj

common_SNPs <- intersect(all_IVs_adj, ld_snps_adj)
cat("Common SNPs between merged data and LD matrix:\n")
print(common_SNPs)
cat("Number of common SNPs:", length(common_SNPs), "\n")
if(length(common_SNPs) == 0){
  stop("No common SNPs found between merged data and LD matrix. Please check your identifiers.")
}
merged_data <- merged_data[SNP %in% common_SNPs]
setorder(merged_data, SNP)
ld_matrix <- ld_matrix_result[merged_data$SNP, merged_data$SNP, drop = FALSE]

# Regularize the LD matrix to avoid singularity by adding a small constant to the diagonal
epsilon <- 1e-5
ld_matrix <- ld_matrix + diag(epsilon, nrow(ld_matrix))

# -------------------------------
# 6. Define Exposure and Outcome Columns for MR Input
# -------------------------------
# In the merged data, we assume the exposure associations are:
# For Expression: "slope_expr" and "slope_se_expr"
# For Splicing:   "slope_splice" and "slope_se_splice"
listBeta <- c("slope_expr", "slope_splice")
listSE   <- c("slope_se_expr", "slope_se_splice")
if (!("beta.outcome" %in% names(merged_data)) || !("se.outcome" %in% names(merged_data))){
  stop("Merged data does not contain outcome columns 'beta.outcome' and 'se.outcome'.")
}
listExposure <- c("Expression", "Splicing")
outcome <- "BMI"
seuil <- 0.05 / length(unique(c(eqtl_data$gene_symbol, sqtl_data$gene_symbol)))

# -------------------------------
# 7. Create the Multivariable MR Input Object
# -------------------------------
# Convert to numeric and check for NAs, zeros, or infinite values
required_cols <- c(listBeta, listSE, "beta.outcome", "se.outcome")
if (!all(required_cols %in% names(merged_data))) {
  stop("Missing required columns in merged_data: ", 
       paste(required_cols[!required_cols %in% names(merged_data)], collapse = ", "))
}
merged_data[, (required_cols) := lapply(.SD, as.numeric), .SDcols = required_cols]

# Check for NA values
if (any(sapply(merged_data[, ..required_cols], function(x) any(is.na(x))))) {
  cat("NA values found in data:\n")
  print(colSums(is.na(merged_data[, ..required_cols])))
}

# Check for zero standard errors
se_cols <- c("slope_se_expr", "slope_se_splice", "se.outcome")
if (any(sapply(merged_data[, ..se_cols], function(x) any(x == 0, na.rm = TRUE)))) {
  cat("Zero standard errors found in columns:\n")
  zero_se_cols <- se_cols[sapply(merged_data[, ..se_cols], function(x) any(x == 0, na.rm = TRUE))]
  print(zero_se_cols)
  print(sapply(merged_data[, ..zero_se_cols], function(x) sum(x == 0, na.rm = TRUE)))
} else {
  cat("No zero standard errors found.\n")
}

# Check for infinite values
if (any(sapply(merged_data[, ..required_cols], function(x) any(is.infinite(x))))) {
  stop("Infinite values found in data.")
}
cat("Summary of exposure and outcome data:\n")
print(summary(merged_data[, ..required_cols]))

# Check LD matrix properties
cat("LD matrix dimensions:", dim(ld_matrix), "\n")
if (any(is.na(ld_matrix))) stop("LD matrix contains NA values.")
if (!isSymmetric(ld_matrix)) stop("LD matrix is not symmetric.")
if (!all(diag(ld_matrix) >= 1 - epsilon & diag(ld_matrix) <= 1 + epsilon)) {
  cat("LD matrix diagonal elements are not approximately 1:\n")
  print(diag(ld_matrix))
}
max_corr <- max(abs(ld_matrix[upper.tri(ld_matrix)]))
cat("Maximum off-diagonal correlation in LD matrix:", max_corr, "\n")
if (max_corr > 0.99) cat("Warning: Very high correlations may cause numerical instability.\n")
eig <- eigen(ld_matrix)$values
if (any(eig <= 0)) {
  cat("LD matrix is not positive definite. Increasing epsilon.\n")
  epsilon <- 1e-4  # Try a larger value
  ld_matrix <- ld_matrix_result[merged_data$SNP, merged_data$SNP, drop = FALSE] + diag(epsilon, nrow(ld_matrix))
}

# Create MR input object
mr_input_object <- mr_mvinput(
  bx = as.matrix(merged_data[, ..listBeta]),
  bxse = as.matrix(merged_data[, ..listSE]),
  by = merged_data[, beta.outcome],
  byse = merged_data[, se.outcome],
  corr = ld_matrix,
  exposure = listExposure,
  outcome = outcome
)

# -------------------------------
# 8. Run Multivariable MR Analyses
# -------------------------------
cat("Running Multivariable MR Analyses...\n")
ivw_result <- try(mr_mvivw(mr_input_object, alpha = seuil), silent = TRUE)
if (inherits(ivw_result, "try-error") || all(is.na(ivw_result@Estimate))) {
  cat("IVW MR analysis failed or produced NA estimates.\n")
  # Test without correlation
  mr_input_uncorr <- mr_mvinput(
    bx = as.matrix(merged_data[, ..listBeta]),
    bxse = as.matrix(merged_data[, ..listSE]),
    by = merged_data[, beta.outcome],
    byse = merged_data[, se.outcome],
    exposure = listExposure,
    outcome = outcome
  )
  ivw_uncorr <- mr_mvivw(mr_input_uncorr, alpha = seuil)
  cat("IVW MR Results (assuming uncorrelated variants):\n")
  print(ivw_uncorr)
} else {
  cat("IVW MR Analysis Results:\n")
  print(ivw_result)
}

egger_result <- try(mr_mvegger(mr_input_object, alpha = seuil), silent = TRUE)
if (inherits(egger_result, "try-error") || all(is.na(egger_result@Estimate))) {
  cat("MR-Egger analysis failed or produced NA estimates.\n")
} else {
  cat("MR-Egger Analysis Results:\n")
  print(egger_result)
}

# Save results if successful
if (!inherits(ivw_result, "try-error") && !all(is.na(ivw_result@Estimate))) {
  write.table(ivw_result, file = "IVW_MR_result.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
}
if (!inherits(egger_result, "try-error") && !all(is.na(egger_result@Estimate))) {
  write.table(egger_result, file = "MR_Egger_result.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
}
cat("MR analysis completed. Check output files if generated.\n")
