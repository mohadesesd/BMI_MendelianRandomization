# -----------------------------------------------------------
# Load Required Libraries
# -----------------------------------------------------------
library(data.table)

# -----------------------------------------------------------
# 1. Read the Processed QTL Files
# -----------------------------------------------------------
# Update the paths if necessary
sqtls <- fread("sQTL_with_SNP.tsv")
eqtls <- fread("eQTL_with_SNP.tsv")

# Inspect the first few rows (optional)
cat("First few rows of sQTL data:\n")
print(head(sqtls))
cat("First few rows of eQTL data:\n")
print(head(eqtls))

# -----------------------------------------------------------
# 2. Identify the Strongest Instrument (Highest Fstat)
# -----------------------------------------------------------
# For the sQTL data:
if("Fstat" %in% names(sqtls)){
  strongest_sqtl <- sqtls[which.max(Fstat)]
  cat("Strongest sQTL IV (by highest Fstat):\n")
  print(strongest_sqtl)
} else {
  cat("Column Fstat not found in sQTL data. Please compute Fstat first.\n")
}

# For the eQTL data:
if("Fstat" %in% names(eqtls)){
  strongest_eqtl <- eqtls[which.max(Fstat)]
  cat("Strongest eQTL IV (by highest Fstat):\n")
  print(strongest_eqtl)
} else {
  cat("Column Fstat not found in eQTL data. Please compute Fstat first.\n")
}

# -----------------------------------------------------------
# 3. (Optional) Save the Results as TSV Files
# -----------------------------------------------------------
fwrite(strongest_sqtl, file = "Strongest_sQTL_IV.tsv", sep = "\t")
fwrite(strongest_eqtl, file = "Strongest_eQTL_IV.tsv", sep = "\t")

cat("The strongest IVs for sQTLs and eQTLs have been saved as 'Strongest_sQTL_IV.tsv' and 'Strongest_eQTL_IV.tsv'.\n")
