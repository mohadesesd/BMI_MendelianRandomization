# -----------------------------------------------------------
# Load Required Libraries
# -----------------------------------------------------------
library(GenomicRanges)   # For handling genomic regions
library(rtracklayer)     # For liftover functions and chain file import
library(data.table)      # For fast data I/O
library(BEDMatrix)       # To read PLINK binary BED files
library(TwoSampleMR)
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(reshape2)
library(LDlinkR)
# -----------------------------------------------------------
# Step 1: Define the Target Gene Region (Already in GRCh38)
# -----------------------------------------------------------
# Define the target gene region (for example, KNOP1) in build38 coordinates.
gene_region_hg38 <- GRanges(seqnames = "chr16", 
                            ranges = IRanges(start = 19701937, end = 19718235))
cat("Target gene region (build38):\n")
print(gene_region_hg38)

# -----------------------------------------------------------
# Step 2: Load the 1000G Data from a Binary BED File
# -----------------------------------------------------------
# Provide the file paths for your binary BED and companion BIM files.
bed_file_path <- "./KNOP1_Adipose_vis/Chr16.EUR_subset.bed"    
bim_file_path <- "./KNOP1_Adipose_vis/Chr16.EUR_subset.bim"      

# Load the BED file as a matrix-like object.
bed_data <- BEDMatrix(bed_file_path)

# Load the companion BIM file, which contains SNP metadata.
bim_data <- fread(bim_file_path, header = FALSE)
setnames(bim_data, old = names(bim_data), 
         new = c("chr", "snp", "gd", "bp", "allele1", "allele2"))
cat("First few lines of BIM data:\n")
print(head(bim_data))

# -----------------------------------------------------------
# Step 3: Liftover the BIM Coordinates from hg19 to hg38
# -----------------------------------------------------------
# Assume that the BIM file SNP positions are in hg19.
# Create a GRanges object from the BIM file.
# (Here we assume that the BIM chromosomes do not include "chr"; we add it.)
bim_gr <- GRanges(seqnames = paste0("chr", bim_data$chr),
                  ranges = IRanges(start = bim_data$bp, end = bim_data$bp))

# Import the chain file (ensure the file "hg19ToHg38.over.chain" is in your working directory or provide a full path).
chain <- import.chain("hg19ToHg38.over.chain")

# Perform liftover on the BIM positions.
lifted_bim <- liftOver(bim_gr, chain)

# Use sapply to create a numeric vector of liftover start positions:
lifted_bim_bp <- sapply(lifted_bim, function(x) {
  if (length(x) > 0) {
    start(x)[1]  # Take the first mapping if multiple exist.
  } else {
    NA_integer_  # Return NA if no mapping is available.
  }
})

# Add the liftover positions to the BIM data as a new column 'bp_hg38'.
bim_data[, bp_hg38 := lifted_bim_bp]
cat("BIM data with liftover coordinates (first few lines):\n")
print(head(bim_data))

# -----------------------------------------------------------
# Step 4: Filter the 1000G Data for the Target Gene Region
# -----------------------------------------------------------
# Define the target region parameters from the gene_region_hg38.
target_chr   <- as.character(seqnames(gene_region_hg38))  # Expected "chr16"
target_start <- start(gene_region_hg38)
target_end   <- end(gene_region_hg38)
cat("Target region (build38):", target_chr, target_start, target_end, "\n")

# Adjust target_chr if the BIM file does not include the "chr" prefix.
if (!any(grepl("^chr", bim_data$chr))) {
  target_chr <- gsub("chr", "", target_chr)  # Remove "chr", e.g., "chr16" becomes "16"
  cat("Adjusted target chromosome:", target_chr, "\n")
}

# Filter the BIM data to include only SNPs with liftover positions that fall within the target region.
bim_filtered <- subset(bim_data, 
                       chr == target_chr & 
                       bp_hg38 >= target_start & 
                       bp_hg38 <= target_end)

if (nrow(bim_filtered) == 0) {
  stop("No SNPs found in the BIM file for the target gene region. Please check the chromosome naming or coordinates.")
} else {
  cat("Number of SNPs in target region:", nrow(bim_filtered), "\n")
  print(head(bim_filtered))
}

# Identify the indices of SNPs that fall in the target region.
snp_indices <- which(bim_data$chr == target_chr & bim_data$bp_hg38 >= target_start & bim_data$bp_hg38 <= target_end)
if (length(snp_indices) == 0) {
  stop("No matching SNP indices found. Verify that the BIM file and the liftover coordinates use matching conventions.")
}

# Subset the genotype matrix from BEDMatrix using these SNP indices.
bed_data_region <- bed_data[, snp_indices]

cat("Dimensions of the filtered genotype matrix:\n")
print(dim(bed_data_region))
cat("Subset of genotype matrix (first 5 rows and 5 SNPs):\n")
print(bed_data_region[1:5, 1:5])


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
sqtl_file_path <- "/lustre03/project/6050805/SHARED_FULL/SHARED/GTEx/sQTL_GTEx_v8/GTEx_Analysis_v8_sQTL_all_associations/All/Adipose/Adipose_Subcutaneous.v8.cis_sqtl.all_pairs.chr16.parquet"

# -----------------------------------------------------------
# Open the Parquet dataset using the lazy API (without memory_map argument)
# -----------------------------------------------------------
ds <- open_dataset(sqtl_file_path, format = "parquet")

# (Optional) Print the dataset schema to see available columns.
print(ds)

available_cols <- names(ds)
cat("Columns to be selected:\n")
print(available_cols)

# -----------------------------------------------------------
# Filter the dataset lazily to only include rows for chromosome 16.
# (We assume that the variant_id field encodes chromosome information.
#  Here we filter rows where the variant_id starts with "chr16".)
ds_filtered <- ds %>% 
  filter(substr(variant_id, 1, 5) == "chr16") %>% 
  select(all_of(available_cols))

# Collect the filtered data into memory as a data.table.
sqtls_data <- as.data.table(ds_filtered %>% collect())
rm(ds, ds_filtered)  # Free up memory from lazy objects

# Preview the first few rows of the loaded data.
cat("First few rows of the filtered sQTL data:\n")
print(head(sqtls_data))

# -----------------------------------------------------------
# Extract Chromosome and Position Information from the variant_id Column
# -----------------------------------------------------------
# The variant_id field is expected to look like: "chr1_64764_C_T_b38"
# We split the string into components: chromosome, position, ref, alt, and build.
sqtls_data[, c("chr", "pos", "ref", "alt", "build") := tstrsplit(variant_id, "_", fixed = TRUE)]
# Convert the 'pos' column to numeric.
sqtls_data[, pos := as.numeric(pos)]

cat("First few rows after extracting 'chr' and 'pos':\n")
print(head(sqtls_data))

cat("Target gene region (build38):\n")
print(gene_region_hg38)

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

#Finding R2 for all values
sqtls_region$R2 <- (get_r_from_bsen(sqtls_region$slope,sqtls_region$slope_se,469))^2 # CHANGE FOR EACH TISSUE

#Finding the Fstat
sqtls_region$Fstat <- (sqtls_region$R2/1)/((1-sqtls_region$R2)/(469-1-1))

#Filtering for Fstat
sqtls_region<- sqtls_region[sqtls_region$Fstat > 10,]

# Split the phenotype_id on the colon delimiter and extract the 5th element
# This stores the full gene id (with version) in a new column called 'gene_full'
sqtls_region[, gene_full := tstrsplit(phenotype_id, ":", fixed = TRUE)[[5]]]

# Remove the version number (dot and following digits)
# This creates a new column 'gene_symbol'
sqtls_region[, gene_symbol := sub("\\..*$", "", gene_full)]

seuil <- 0.05/length(unique(sqtls_region$gene_symbol))

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

#Finding R2 for all values
eqtls_region$R2 <- (get_r_from_bsen(eqtls_region$slope,eqtls_region$slope_se,469))^2 # CHANGE FOR EACH TISSUE

#Finding the Fstat
eqtls_region$Fstat <- (eqtls_region$R2/1)/((1-eqtls_region$R2)/(469-1-1))

#Filtering for Fstat
eqtls_region<- eqtls_region[eqtls_region$Fstat > 10,]


# Remove the version number (dot and following digits)
# This creates a new column 'gene_symbol'
eqtls_region[, gene_symbol := sub("\\..*$", "", gene_id)]

seuil <- 0.05/length(unique(eqtls_region$gene_symbol))

eqtls_region <- eqtls_region[eqtls_region$pval_nominal < seuil, ]

cat("Filtered eQTL data in the  KNOP1 region (build38):\n")
print(head(eqtls_region))

# -----------------------------------------------------------

library(data.table)

# 1. Save the Filtered QTL Data as TSV Files

# Save the filtered sQTL data (sqtls_region)
fwrite(sqtls_region, file = "Filtered_sQTLs_KNOP1.tsv", sep = "\t")
cat("Filtered sQTLs saved as 'Filtered_sQTLs_KNOP1.tsv'\n")

# Save the filtered eQTL data (eqtls_region)
fwrite(eqtls_region, file = "Filtered_eQTLs_KNOP1.tsv", sep = "\t")
cat("Filtered eQTLs saved as 'Filtered_eQTLs_KNOP1.tsv'\n")


# 2. Save the Subsetted (Lifted‐over) BIM/Bed Files in PLINK Format

# We already have a filtered BIM object (bim_filtered) that contains the SNPs in the target region.
# First, extract the SNP IDs from bim_filtered.
snp_list <- bim_filtered$snp

# Write the SNP IDs to a file (one SNP per line, no header).
# (The file "snp_list.txt" will be used by PLINK to extract these SNPs.)
fwrite(data.table(snp = snp_list), file = "snp_list.txt", sep = "\t", col.names = FALSE)
cat("SNP list saved as 'snp_list.txt'\n")

# Next, call PLINK from within R to subset the original binary dataset.
# (Assuming your original binary files are named "Chr16.EUR_subset.bed", "Chr16.EUR_subset.bim", "Chr16.EUR_subset.fam".)
# The following PLINK command extracts the SNPs in the list and writes a new binary dataset.
plink_command <- "plink --bfile Chr16.EUR_subset --extract snp_list.txt --make-bed --out Subsetted_Lifted_KNOP1"
cat("Running PLINK command:\n", plink_command, "\n")
system(plink_command)

cat("Subsetted PLINK files (BED, BIM, and FAM) have been saved with the prefix 'Subsetted_Lifted_KNOP1'.\n")
