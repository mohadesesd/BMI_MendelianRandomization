# -----------------------------------------------------------
# IV_selection_and_LD_pruning_separately.R
#
# This script processes exposure (QTL) data separately for:
#  - Expression (eQTL) data, and
#  - Splicing (sQTL) data.
#
# It reads the QTL files, cleans and deduplicates them,
# extracts instrument strength values (Fstat), loads a provided LD matrix,
# harmonizes SNP IDs, performs LD pruning (using an r² threshold),
# and then saves two separate files with the final instrument lists.
# -----------------------------------------------------------

# Load required libraries
library(data.table)  # For fast data manipulation

# -------------------------------
# 1. Read and Process Expression (eQTL) Data
# -------------------------------
eqtl_data <- fread("eQTL_with_SNP.tsv")
# Remove rows with missing SNPs
eqtl_data <- eqtl_data[!is.na(SNP) & SNP != ""]
# Deduplicate by SNP (if there are duplicates, keep the first)
eqtl_data <- unique(eqtl_data, by = "SNP")
# (Assumption: The eQTL file has a column "Fstat" representing the instrument strength.)
expr_iv <- eqtl_data[, .(SNP, Fstat)]
# Ensure SNP identifiers have the "rs" prefix
expr_iv[, SNP := ifelse(grepl("^rs", SNP), SNP, paste0("rs", SNP))]
cat("Number of unique eQTL instruments:", nrow(expr_iv), "\n")
cat("eQTL SNPs:\n")
print(expr_iv$SNP)

# -------------------------------
# 2. Read and Process Splicing (sQTL) Data
# -------------------------------
sqtl_data <- fread("sQTL_with_SNP.tsv")
sqtl_data <- sqtl_data[!is.na(SNP) & SNP != ""]
sqtl_data <- unique(sqtl_data, by = "SNP")
# (Assumption: sQTL file also has a column "Fstat" for instrument strength.)
splice_iv <- sqtl_data[, .(SNP, Fstat)]
splice_iv[, SNP := ifelse(grepl("^rs", SNP), SNP, paste0("rs", SNP))]
cat("Number of unique sQTL instruments:", nrow(splice_iv), "\n")
cat("sQTL SNPs:\n")
print(splice_iv$SNP)

# -------------------------------
# 3. Load and Process the LD Matrix
# -------------------------------
# Read LD matrix from file.
ld_matrix <- read.table("LD_matrix_result.tsv",
                        header = TRUE,
                        sep = "\t",
                        row.names = 1,
                        stringsAsFactors = FALSE,
                        fill = TRUE,
                        check.names = FALSE)
ld_matrix <- as.matrix(ld_matrix)
cat("Original LD matrix dimensions:", dim(ld_matrix), "\n")
# Ensure LD matrix is square: if not, take the first n_rows columns.
n_rows <- nrow(ld_matrix)
n_cols <- ncol(ld_matrix)
if(n_rows != n_cols){
  cat("LD matrix is not square. Subsetting to the first", n_rows, "columns.\n")
  ld_matrix <- ld_matrix[, 1:n_rows]
  cat("New LD matrix dimensions:", dim(ld_matrix), "\n")
}
# Harmonize LD matrix SNP IDs
ld_snps <- rownames(ld_matrix)
ld_snps_adj <- ifelse(grepl("^rs", ld_snps), ld_snps, paste0("rs", ld_snps))
rownames(ld_matrix) <- ld_snps_adj
colnames(ld_matrix) <- ld_snps_adj
cat("Adjusted SNP IDs from LD matrix:\n")
print(ld_snps_adj)

# -------------------------------
# 4. Subset LD Matrix to Overlap with Instruments for Each Exposure
# -------------------------------

# For Expression instruments:
common_expr <- intersect(expr_iv$SNP, ld_snps_adj)
expr_iv_final <- expr_iv[SNP %in% common_expr]
cat("Number of eQTL instruments overlapping with LD matrix:", nrow(expr_iv_final), "\n")
ld_expr <- ld_matrix[common_expr, common_expr, drop = FALSE]

# For Splicing instruments:
common_splice <- intersect(splice_iv$SNP, ld_snps_adj)
splice_iv_final <- splice_iv[SNP %in% common_splice]
cat("Number of sQTL instruments overlapping with LD matrix:", nrow(splice_iv_final), "\n")
ld_splice <- ld_matrix[common_splice, common_splice, drop = FALSE]

# -------------------------------
# 5. LD-Based Pruning Function
# -------------------------------
ld_prune <- function(iv_dt, ld_mat, threshold = 0.9) {
  snp_ids <- iv_dt$SNP
  n <- length(snp_ids)
  keep <- rep(TRUE, n)
  names(keep) <- snp_ids
  for(i in 1:(n - 1)){
    if(!keep[i]) next  # Skip if SNP i already pruned
    for(j in (i + 1):n){
      if(!keep[j]) next  # Skip if SNP j already pruned
      if(abs(ld_mat[i, j]) > threshold){
        # Compare F-statistics
        F_i <- iv_dt[SNP == snp_ids[i], Fstat]
        F_j <- iv_dt[SNP == snp_ids[j], Fstat]
        # Use -Inf if missing
        if(is.na(F_i)) F_i <- -Inf
        if(is.na(F_j)) F_j <- -Inf
        if(F_i >= F_j){
          keep[j] <- FALSE
        } else {
          keep[i] <- FALSE
        }
      }
    }
  }
  final_snps <- snp_ids[keep]
  return(final_snps)
}

# -------------------------------
# 6. Prune Instruments Separately for Expression and Splicing
# -------------------------------
ld_threshold <- 0.9

final_expr_snps <- ld_prune(expr_iv_final, ld_expr, threshold = ld_threshold)
cat("Final retained eQTL (Expression) SNPs:\n")
print(final_expr_snps)

final_splice_snps <- ld_prune(splice_iv_final, ld_splice, threshold = ld_threshold)
cat("Final retained sQTL (Splicing) SNPs:\n")
print(final_splice_snps)

# -------------------------------
# 7. Save the Final Instrument Files
# -------------------------------
fwrite(expr_iv_final[SNP %in% final_expr_snps], file = "Final_IVs_expr.tsv", sep = "\t")
fwrite(splice_iv_final[SNP %in% final_splice_snps], file = "Final_IVs_splice.tsv", sep = "\t")
cat("Final IV data saved as 'Final_IVs_expr.tsv' and 'Final_IVs_splice.tsv'.\n")
# Load required library
library(data.table)

# Read the final IVs list (assuming it contains a column named "SNP")
#final_iv <- fread("Final_IVs.tsv")
# Make sure the SNP identifiers are stored in a vector
#final_SNPs <- final_iv$SNP
#cat("Number of final IV SNPs:", length(final_SNPs), "\n")

# ----------------------------
# Subset the eQTL data
# ----------------------------

# Read the eQTL data file
eqtl_data <- fread("eQTL_with_SNP.tsv")
# Optionally, inspect the data
cat("Original eQTL data dimensions:", dim(eqtl_data), "\n")

# Subset the eQTL data to keep only rows where SNP is in the final IV list
eqtl_final <- eqtl_data[SNP %in% final_expr_snps]
cat("Final eQTL IVs (subset) dimensions:", dim(eqtl_final), "\n")

# Save the subsetted eQTL data
fwrite(eqtl_final, file = "Final_IVs_eQTL.tsv", sep = "\t")
cat("Final eQTL IV data saved as 'Final_IVs_eQTL.tsv'.\n")

# ----------------------------
# Subset the sQTL data
# ----------------------------

# Read the sQTL data file
sqtl_data <- fread("sQTL_with_SNP.tsv")
cat("Original sQTL data dimensions:", dim(sqtl_data), "\n")

# Subset the sQTL data to keep only rows where SNP is in the final IV list
sqtl_final <- sqtl_data[SNP %in% final_splice_snps]
cat("Final sQTL IVs (subset) dimensions:", dim(sqtl_final), "\n")

# Save the subsetted sQTL data
fwrite(sqtl_final, file = "Final_IVs_sQTL.tsv", sep = "\t")
cat("Final sQTL IV data saved as 'Final_IVs_sQTL.tsv'.\n")
