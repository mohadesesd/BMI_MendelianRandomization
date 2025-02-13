
# Header ---------------------------------------------------------------------

# AUTHOR: Mohadese Sayahian Dehkordi
# PURPOSE: Liftover outcome from hg19 to 38
# NOTES: 
#   *Inputs:   - Outcome data 
#
#   *Outputs:  - liftover data

# Load required libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("rtracklayer")
library(rtracklayer)
library(GenomicRanges)
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(reshape2)
# Import the chain file (adjust the path to your chain file)
chain <- import.chain("/lustre03/project/6050805/sayahian/hg19ToHg38.over.chain")

# Create a GRanges object from coordinates
freadf <- function(df){return(as.data.frame(fread(df)))}

# Load dataset ----------------------------------------------------------------
df <- freadf("/lustre03/project/6050805/SHARED_FULL/SHARED/GWAS/BMI/pediatric_BMI/Vogelezang_2020/GCST90002409_buildGRCh37.tsv")
print(head(df))
## Ensure chromosome names have the "chr" prefix
df$chromosome <- ifelse(grepl("^chr", df$chromosome),
                        df$chromosome,
                        paste0("chr", df$chromosome))
## Ensure base_pair locations are numeric
df$base_pair_location <- as.numeric(as.character(df$base_pair_location))
### Check the number of rows in each column
cat("Length of chromosome:", length(df$chromosome), "\n")
cat("Length of base_pair_location:", length(df$base_pair_location), "\n")
# Create GRanges object from the data frame.
## Since variants are single nucleotide positions, we use the same value for start and end.
gr <- GRanges(
  seqnames = df$chromosome,
  ranges = IRanges(start = df$base_pair_location, end = df$base_pair_location),
  variant_id = df$variant_id,
  chromosome = df$chromosome, 
  base_pair_location = df$base_pair_location, 
  effect_allele = df$effect_allele, 
  other_allele = df$other_allele, 
  beta = df$beta, 
  standard_error = df$standard_error,
  p_value = df$p_value, 
  TotalSampleSize = df$TotalSampleSize
)
print(head(gr))
# Perform liftover
lifted <- liftOver(gr, chain)

# Unlist if expecting one-to-one mappings
lifted_gr <- unlist(lifted)

# Convert results to a data frame and view
lifted_df <- as.data.frame(lifted_gr)
# Orgnize columns
lifted_df <- cbind(lifted_df[,c(6,7,2,9,10,11,12,13,14)])
colnames(lifted_df) <- c("variant_id", "chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", "standard_error", "p_value", "TotalSampleSize")

print(head(lifted_df))

write.table(lifted_df, file = "GCST90002409_buildGRCh38.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
