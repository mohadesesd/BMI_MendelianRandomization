# -----------------------------------------------------------
# Load Required Libraries
# -----------------------------------------------------------
library(LDlinkR)
library(data.table)

# -----------------------------------------------------------
# 1. Read the eQTL and sQTL TSV Files
# -----------------------------------------------------------
# Update the file paths to where your TSV files are located.
eqtl_data <- fread("eQTL_with_SNP.tsv")
sqtl_data <- fread("sQTL_with_SNP.tsv")

# -----------------------------------------------------------
# 2. Prepare the List of SNPs (IVs)
# -----------------------------------------------------------
# Assuming that both files contain a column "SNP" with the rsIDs,
# combine them and keep only unique non-missing values.
all_SNPs <- unique(c(eqtl_data$SNP, sqtl_data$SNP))
# Remove NA or empty strings
all_SNPs <- all_SNPs[!is.na(all_SNPs) & all_SNPs != ""]

cat("Total number of unique SNPs (IVs):", length(all_SNPs), "\n")

# -----------------------------------------------------------
# 3. Set Up LDlinkR Parameters and Compute the LD Matrix
# -----------------------------------------------------------
# Replace with your valid LDlink API token.
token <- "65d9498c4838"

# Compute the LD matrix for the provided SNP list.
# The pop parameter is set to "EUR" for European population.
ld_matrix_result <- LDmatrix(
  snps = all_SNPs,
  pop = "EUR",         # European population (adjust if needed)
  r2d = "r2",          # Return r² values (change to "dprime" for D')
  token = token
)

# -----------------------------------------------------------
# 4. Save the LD Matrix as a TSV File
# -----------------------------------------------------------
# Save ld_matrix_result as a TSV file (tab-separated values)
write.table(ld_matrix_result, file = "LD_matrix_result.tsv", 
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

cat("LD matrix saved as 'LD_matrix_result.tsv'\n")
