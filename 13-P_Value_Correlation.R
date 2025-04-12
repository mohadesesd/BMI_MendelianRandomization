#!/usr/bin/env Rscript
library(dplyr)
library(readr)
library(stringr)

#-------------------------------------------
# Set directories for cleaned MR files
#-------------------------------------------
adult_dir     <- "./LCL_MR/Cleaned_MR2/"      # Directory containing adult MR files (e.g., cleaned_MR_chr*.rds)
pediatric_dir <- "./LCL_MR/Cleaned_Raw_MR/"   # Directory containing pediatric MR files
output_csv    <- "./LCL_MR/correlation_results.csv" # File to save correlation results

#-------------------------------------------
# Step 1: Combine all adult MR data and extract sQTL END from id.exposure
#-------------------------------------------
adult_files <- list.files(path = adult_dir, pattern = "^cleaned_MR_chr.*\\.rds$", full.names = TRUE)
if (length(adult_files) == 0) stop("No adult MR files found in the directory.")

adult_list <- lapply(adult_files, readRDS)
adult_data <- bind_rows(adult_list)
cat("Adult MR data combined:", nrow(adult_data), "rows.\n")

# Extract sQTL end position from id.exposure
# Assumes id.exposure is in the format "chrX:start:end:clu_xxxx:gene_id.version"
adult_data <- adult_data %>%
  mutate(sQTL_end = as.numeric(str_split(id.exposure, ":", simplify = TRUE)[,3]))

#-------------------------------------------
# Step 2: Combine all pediatric MR data and extract sQTL END
#-------------------------------------------
ped_files <- list.files(path = pediatric_dir, pattern = "^cleaned_MR_chr.*\\.rds$", full.names = TRUE)
if (length(ped_files) == 0) stop("No pediatric MR files found in the directory.")

ped_list <- lapply(ped_files, readRDS)
ped_data <- bind_rows(ped_list)
cat("Pediatric MR data combined:", nrow(ped_data), "rows.\n")

# Extract sQTL end for pediatric data
ped_data <- ped_data %>%
  mutate(sQTL_end = as.numeric(str_split(id.exposure, ":", simplify = TRUE)[,3]))

#-------------------------------------------
# Step 3: Aggregate data by gene_id and chr_final
#-------------------------------------------
# Aggregate to one record per gene by taking the median p-value and median sQTL end.
adult_data_agg <- adult_data %>%
  group_by(gene_id, chr_final) %>%
  summarize(pval.adult = median(pval, na.rm = TRUE),
            sQTL_end.adult = median(sQTL_end, na.rm = TRUE),
            n_adult = n(),
            .groups = "drop")

ped_data_agg <- ped_data %>%
  group_by(gene_id, chr_final) %>%
  summarize(pval.ped = median(pval, na.rm = TRUE),
            sQTL_end.ped = median(sQTL_end, na.rm = TRUE),
            n_ped = n(),
            .groups = "drop")

cat("Adult data aggregated to", nrow(adult_data_agg), "unique gene entries.\n")
cat("Pediatric data aggregated to", nrow(ped_data_agg), "unique gene entries.\n")

#-------------------------------------------
# Step 4: Merge aggregated data by gene_id and chr_final
#-------------------------------------------
merged_data <- inner_join(adult_data_agg, ped_data_agg, by = c("gene_id", "chr_final"))
cat("Merged data contains:", nrow(merged_data), "rows.\n")

#-------------------------------------------
# Step 5: Exclude records in the MHC region (chr6: 28,510,120 to 33,480,577) based on sQTL END
#-------------------------------------------
merged_data <- merged_data %>%
  filter(!(chr_final == "chr6" & sQTL_end.adult >= 28510120 & sQTL_end.adult <= 33480577))

cat("After excluding MHC region, merged data contains:", nrow(merged_data), "rows.\n")

#-------------------------------------------
# Step 6: Compute Pearson correlation between aggregated p-values
#-------------------------------------------
if (!("pval.adult" %in% colnames(merged_data)) || !("pval.ped" %in% colnames(merged_data))) {
  stop("p-value columns not found in the merged dataset.")
}

cor_val <- cor(merged_data$pval.adult, merged_data$pval.ped, method = "pearson", use = "complete.obs")
cat("Pearson correlation between adult and pediatric aggregated p-values:", cor_val, "\n")

#-------------------------------------------
# Step 7: Save the correlation results to CSV
#-------------------------------------------
cor_results <- data.frame(
  overall_correlation = cor_val,
  n_records = nrow(merged_data)
)

write_csv(cor_results, output_csv)
cat("Correlation results saved to", output_csv, "\n")
