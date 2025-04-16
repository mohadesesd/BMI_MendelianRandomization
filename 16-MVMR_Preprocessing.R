# -----------------------------------------------------------
# Load Required Libraries
# -----------------------------------------------------------
library(GenomicRanges)   # For handling genomic regions
library(rtracklayer)     # For liftover functions and chain file import
library(data.table)      # For fast data I/O
library(BEDMatrix)       # To read PLINK binary BED files
library(TwoSampleMR)
library(dplyr)
library(tidyr)
library(stringr)
library(reshape2)
library(LDlinkR)

# -----------------------------------------------------------
# Step 1: Define the Target Gene Region (Already in GRCh38)
# -----------------------------------------------------------
# Define the target gene region (for example, KNOP1) in build38 coordinates.
gene_region_hg38 <- GRanges(seqnames = "chr2", 
                            ranges = IRanges(start = 24819169, end = 24920237))
cat("Target gene region (build38):\n")
print(gene_region_hg38)

# -----------------------------------------------------------
# Step 5: Extract sQTL Data from GTEx for the Gene Region
# -----------------------------------------------------------
# Step: Load and Filter GTEx sQTL Data from a gzipped TXT File (No Liftover Needed)
# -----------------------------------------------------------
if (!requireNamespace("R.utils", quietly = TRUE)) {
  install.packages("R.utils", repos = "https://cloud.r-project.org")
}

# -----------------------------------------------------------
# Load Required Libraries
# -----------------------------------------------------------
library(arrow)         # For reading Parquet files via the dataset API
library(dplyr)         # For lazy filtering/selection
library(data.table)    # For fast in-memory table operations
library(GenomicRanges) # For defining genomic regions

# -----------------------------------------------------------
# Specify the file path to your GTEx sQTL Parquet file for Visceral Adipose.
# -----------------------------------------------------------
sqtl_file_path <- "/lustre03/project/6050805/SHARED_FULL/SHARED/GTEx/sQTL_GTEx_v8/GTEx_Analysis_v8_sQTL_all_associations/All/Adipose/Adipose_Subcutaneous.v8.cis_sqtl.all_pairs.chr2.parquet"

# -----------------------------------------------------------
# Open the Parquet dataset using the lazy API (without memory_map argument)
# -----------------------------------------------------------
ds <- open_dataset(sqtl_file_path, format = "parquet")

# (Optional) Print the dataset schema to see available columns.
sqtls_data <- as.data.frame(ds)

available_cols <- names(ds)
cat("Columns to be selected:\n")
print(available_cols)

# -----------------------------------------------------------
# Filter the dataset lazily to only include rows for chromosome 16.
# (We assume that the variant_id field encodes chromosome information.
#  Here we filter rows where the variant_id starts with "chr16".)
# Remove leading/trailing spaces and then filter
sqtls_data <- as.data.table(sqtls_data) # Convert to data.table
sqtls_data[, c("chr", "pos", "ref", "alt", "build") := tstrsplit(variant_id, "_", fixed = TRUE)]
sqtls_data[, pos := as.numeric(pos)] # Convert the 'pos' column to numeric

cat("First few rows after extracting 'chr' and 'pos':\n")
print(head(sqtls_data))

# -----------------------------------------------------------
# Filtering sQTL Data for the KNOP1 Region
# -----------------------------------------------------------
target_chr   <- as.character(seqnames(gene_region_hg38))
target_start <- start(gene_region_hg38)
target_end   <- end(gene_region_hg38)
cat("Target region:", target_chr, target_start, target_end, "\n")

# If needed, adjust target_chr based on how your sQTL data stores the chromosome.
# For example, if your data include "chr16" then it's fine; if not, remove "chr" prefix.
if (!any(grepl("^chr", sqtls_data$chr))) {
  target_chr <- gsub("chr", "", target_chr)
}
cat("Using target chromosome for filtering:", target_chr, "\n")

# Filter the sQTL data based on the extracted 'chr' and 'pos' columns.
sqtls_region <- subset(sqtls_data,
                       chr == target_chr &
                       pos >= target_start &
                       pos <= target_end)

# Finding R2 for all values
sqtls_region$R2 <- (get_r_from_bsen(sqtls_region$slope, sqtls_region$slope_se, 581))^2 # CHANGE FOR EACH TISSUE

# Finding the Fstat
sqtls_region$Fstat <- (sqtls_region$R2 / 1) / ((1 - sqtls_region$R2) / (581 - 1 - 1))

# Filtering for Fstat
sqtls_region <- sqtls_region[sqtls_region$Fstat > 10,]

# Split the phenotype_id on the colon delimiter and extract the 5th element
# This stores the full gene id (with version) in a new column called 'gene_full'
sqtls_region[, gene_full := tstrsplit(phenotype_id, ":", fixed = TRUE)[[5]]]

# Remove the version number (dot and following digits)
# This creates a new column 'gene_symbol'
sqtls_region[, gene_symbol := sub("\\..*$", "", gene_full)]

seuil <- 0.05 / length(unique(sqtls_region$gene_symbol))

sqtls_region <- sqtls_region[sqtls_region$pval_nominal < seuil, ]

cat("Filtered sQTL data in the KNOP1 region (build38):\n")
print(head(sqtls_region))

# -----------------------------------------------------------
# Specify the file path to your GTEx eQTL file for Visceral Adipose.
# -----------------------------------------------------------
eqtl_file_path <- "/lustre03/project/6050805/SHARED_FULL/SHARED/GTEx/eQTL_GTEx_v8/GTEx_Analysis_v8_eQTL_all_associations/All/Adipose/Adipose_Subcutaneous.allpairs.txt.gz"

# -----------------------------------------------------------
# Load the gzipped eQTL data using fread()
# -----------------------------------------------------------
eqtl_data <- fread(eqtl_file_path)

# Preview the first few rows to check the column names.
cat("First few rows of the Visceral Adipose eQTL data:\n")
print(head(eqtl_data))

# -----------------------------------------------------------
# Extract Chromosome and Position Information from the variant_id Column
# -----------------------------------------------------------
# The variant_id field is expected to look like: "chr1_64764_C_T_b38"
# We use tstrsplit() to split the string into components:
# chromosome, position, ref, alt, and build.
eqtl_data[, c("chr", "pos", "ref", "alt", "build") := tstrsplit(variant_id, "_", fixed = TRUE)]
# Convert the 'pos' column to numeric.
eqtl_data[, pos := as.numeric(pos)]

cat("First few rows after extracting 'chr' and 'pos':\n")
print(head(eqtl_data))

cat("Target gene region (build38):\n")
print(gene_region_hg38)

# -----------------------------------------------------------
# Filtering eQTL Data for the KNOP1 Region
# -----------------------------------------------------------
target_chr   <- as.character(seqnames(gene_region_hg38))
target_start <- start(gene_region_hg38)
target_end   <- end(gene_region_hg38)
cat("Target region:", target_chr, target_start, target_end, "\n")

# If eQTL data stores chromosome values with "chr" (e.g., "chr16"), then all is fine.
# Otherwise, if the data do not include the "chr" prefix, adjust target_chr by removing "chr":
if (!any(grepl("^chr", eqtl_data$chr))) {
    target_chr <- gsub("chr", "", target_chr)
}
cat("Using target chromosome for filtering:", target_chr, "\n")

# Filter the eQTL data based on the extracted 'chr' and numeric 'pos'.
eqtls_region <- subset(eqtl_data,
                      chr == target_chr &
                      pos >= target_start &
                      pos <= target_end)

# Finding R2 for all values
eqtls_region$R2 <- (get_r_from_bsen(eqtls_region$slope, eqtls_region$slope_se, 581))^2 # CHANGE FOR EACH TISSUE

# Finding the Fstat
eqtls_region$Fstat <- (eqtls_region$R2 / 1) / ((1 - eqtls_region$R2) / (581 - 1 - 1))

# Filtering for Fstat
eqtls_region <- eqtls_region[eqtls_region$Fstat > 10,]

# Remove the version number (dot and following digits)
# This creates a new column 'gene_symbol'
eqtls_region[, gene_symbol := sub("\\..*$", "", gene_id)]

seuil <- 0.05 / length(unique(eqtls_region$gene_symbol))

eqtls_region <- eqtls_region[eqtls_region$pval_nominal < seuil, ]

cat("Filtered eQTL data in the FKBP11 region (build38):\n")
print(head(eqtls_region))

# -----------------------------------------------------------

library(data.table)

# 1. Save the Filtered QTL Data as TSV Files

# Save the filtered sQTL data (sqtls_region)
fwrite(sqtls_region, file = "Filtered_sQTLs_FKBP11.tsv", sep = "\t")
cat("Filtered sQTLs saved as 'Filtered_sQTLs_FKBP11.tsv'\n")

# Save the filtered eQTL data (eqtls_region)
fwrite(eqtls_region, file = "Filtered_eQTLs_FKBP11.tsv", sep = "\t")
cat("Filtered eQTLs saved as 'Filtered_eQTLs_FKBP11.tsv'\n")
