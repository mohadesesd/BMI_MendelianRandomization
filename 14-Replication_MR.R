# ---------------------------------------------------------------------
# AUTHOR: John Eric Hamilton (Modified)
# PURPOSE: Replication eQTLgen dataset with error handling for empty/messy data,
#          ensuring required columns (including variant_id) are present,
#          and reordering columns so that variant_id and SNP are at the end.
# ---------------------------------------------------------------------

# Load libraries and install packages if needed
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = "https://cloud.r-project.org/")
}
remotes::install_github("MRCIEU/genetics.binaRies")
remotes::install_github("ropensci/rsnps")
library(rsnps)
library(TwoSampleMR)
library(dplyr)

# Loading in the data and Preprocessing ----------------------------------------
df_new <- read.delim("/home/sayahian/projects/def-dman2020/SHARED_FULL/SHARED/GWAS/eQTLgen/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz", 
                     fileEncoding = "UTF-8")
allGenes <- read.delim("allGenes.txt", header = FALSE, fileEncoding = "UTF-8")

eQTL <- df_new[df_new$GeneSymbol %in% allGenes$V1, ]

# Write genes not found in expression data
write.table(allGenes[!allGenes$V1 %in% eQTL$GeneSymbol, ], 
            'genesToProxyNewEQTL.txt', row.names = FALSE, quote = FALSE, col.names = FALSE)

# Calculate beta and SE from Zscore and Pvalue
eQTL$Beta <- eQTL$Zscore / abs(qnorm(eQTL$Pvalue / 2))
eQTL$SE <- eQTL$Beta / eQTL$Zscore

# MAF search -------------------------------------------------------------------
SNP_done <- data.frame()
allSNPs <- eQTL$SNP
outcomeToMerge3 <- eQTL
toSkip <- data.frame()

for (i in 1:nrow(outcomeToMerge3)) {
  tryCatch({
    t <- as.data.frame(ncbi_snp_query(as.character(outcomeToMerge3[i, 'SNP'])))$maf_population[[1]]
    t <- t[t$study == 'dbGaP_PopFreq', ]
    
    t2 <- outcomeToMerge3[i, ]
    t3 <- t[(t$ref_seq == t2$AssessedAllele & t$Minor == t2$OtherAllele) | 
              (t$Minor == t2$AssessedAllele & t$ref_seq == t2$OtherAllele), ]
    
    if (t3$Minor == t2$AssessedAllele) {
      outcomeToMerge3[i, 'EAF'] <- t3$MAF
    } else {
      outcomeToMerge3[i, 'EAF'] <- 1 - t3$MAF
    }
    print(t3$MAF)
    print(i)
  }, error = function(e) {
    toSkip <- rbind(toSkip, outcomeToMerge3[i, 'SNP'])
  })
}

if (!dir.exists('eQTLgen/LCL/')) {
  dir.create('eQTLgen/LCL/', recursive = TRUE)
}
saveRDS(outcomeToMerge3, 'eQTLgen/LCL/eQTL_searchedMAF.rds')
outcomeToMerge_updated <- readRDS('eQTLgen/LCL/eQTL_searchedMAF.rds')

# Clumping (LD) ----------------------------------------------------------------
# Format the exposure data using format_data from TwoSampleMR
formexposure01 <- format_data(
  dat = outcomeToMerge_updated,
  type = "exposure",
  snps = NULL,       # leave as is
  header = TRUE,     # leave as is
  phenotype_col = "GeneSymbol",
  snp_col = "SNP",
  beta_col = "Beta",
  se_col = "SE",
  eaf_col = "EAF", 
  effect_allele_col = "AssessedAllele",
  other_allele_col = "OtherAllele",
  pval_col = "Pvalue",
  units_col = "units",
  ncase_col = "ncase",
  ncontrol_col = "ncontrol",
  samplesize_col = "NrSamples",
  gene_col = "GeneSymbol", 
  id_col = "Gene", 
  min_pval = 1e-200,
  z_col = "z",
  info_col = "info",
  chr_col = "SNPChr",
  pos_col = "SNPPos",
  log_pval = FALSE
)

# Set API option to remove recurring error messages
options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')

# Initialize list to store clumped results for each chromosome
clumped_list <- list()
for (chr in 1:22) {
  dat_chr <- formexposure01[formexposure01$chr.exposure == chr, ]
  
  if (nrow(dat_chr) == 0) {
    message("No data to clump for chromosome ", chr)
    clumped_list[[chr]] <- dat_chr[0, ]
    next
  }
  
  bfile_chr <- paste0("/home/sayahian/projects/def-dman2020/SHARED_FULL/SHARED/1000Genomes/wgs/plink/EUR/Chr", chr, ".EUR")
  
  if (!file.exists(paste0(bfile_chr, ".bed"))) {
    message("PLINK file not found for chromosome ", chr, ": ", bfile_chr, ".bed")
    clumped_list[[chr]] <- dat_chr[0, ]
    next
  }
  
  clumped_chr <- clump_data(
    dat = dat_chr,
    clump_kb = 10000,
    clump_r2 = 0.001,
    clump_p1 = 1,
    clump_p2 = 1,
    pop = "EUR",
    bfile = bfile_chr,
    plink_bin = genetics.binaRies::get_plink_binary()
  )
  clumped_list[[chr]] <- clumped_chr
}

# Combine clumped results into one data frame
clumped_exposure01 <- do.call(rbind, clumped_list)

if (nrow(clumped_exposure01) > 0) {
  # Set the exposure column as id.exposure
  clumped_exposure01$exposure <- clumped_exposure01$id.exposure
  print("hello")
  print(head(clumped_exposure01))
  
} else {
  stop("No clumped data available for any chromosome. Check PLINK file paths and clumping process.")
}

# MR ---------------------------------------------------------------------------
for (i in 1:22) {
  # Extract exposure data for the current chromosome
  eQTL_WB <- clumped_exposure01 %>% filter(chr.exposure == i)
  
  # Check if the data is empty or has only NA/0 values (or if it is a matrix with messy column names)
  if (is.matrix(eQTL_WB) || nrow(eQTL_WB) == 0 || all(is.na(eQTL_WB$SNP))) {
    message("Skipping chromosome ", i, " due to empty or invalid exposure data (matrix/NA).")
    next
  }
  
  # If the data is a matrix, convert it to a data frame and attempt to restore headers.
  # Here we expect the main columns to be in the following order:
  # id.exposure, exposure, beta.exposure, se.exposure, effect_allele.exposure,
  # other_allele.exposure, chr.exposure, pos.exposure, pval.exposure, eaf.exposure,
  # samplesize.exposure, mr_keep.exposure, pval_origin.exposure, variant_id, SNP
  if (is.matrix(eQTL_WB)) {
    eQTL_WB <- as.data.frame(eQTL_WB)
    expected_cols <- c("id.exposure", "exposure", "beta.exposure", "se.exposure", 
                       "effect_allele.exposure", "other_allele.exposure", "chr.exposure", 
                       "pos.exposure", "pval.exposure", "eaf.exposure", "samplesize.exposure", 
                       "mr_keep.exposure", "pval_origin.exposure", "variant_id", "SNP")
    if(ncol(eQTL_WB) == length(expected_cols)){
      colnames(eQTL_WB) <- expected_cols
    } else {
      message("Column number mismatch; attempting to merge with original formatted data for chromosome ", i)
      eQTL_WB <- merge(eQTL_WB, formexposure01, by = "SNP", all.x = TRUE)
    }
  }
  
  # Filter out rows with missing required columns, including variant_id and SNP
  required_cols <- c("id.exposure", "exposure", "beta.exposure", "se.exposure", 
                     "effect_allele.exposure", "other_allele.exposure", "SNP")
  eQTL_WB <- eQTL_WB[complete.cases(eQTL_WB[, required_cols]), ]
  
  if (nrow(eQTL_WB) == 0) {
    message("Skipping chromosome ", i, " after filtering due to no valid exposure data.")
    next
  }
  
  print(head(eQTL_WB))
  
  # Save preprocessed exposure data
  if (!dir.exists('eQTLgen/MAF/Data_MR/LCL/Preprocessed_Data/')) {
    dir.create('eQTLgen/MAF/Data_MR/LCL/Preprocessed_Data/', recursive = TRUE)
  }
  saveRDS(eQTL_WB, sprintf('eQTLgen/MAF/Data_MR/LCL/Preprocessed_Data/eQTL_WB_chr%s.rds', i))
  print(head(eQTL_WB))
  # Load the corresponding outcome data
  GCST4023_Reduced_Renamed <- readRDS(sprintf('LCL_MR/GCST2409/GCST2409_chr%s.rds', i))
  GCST4023_Reduced_Renamed$SNP <- GCST4023_Reduced_Renamed$variant_id
  print(head(GCST4023_Reduced_Renamed))
  if (!"eaf.outcome" %in% names(GCST4023_Reduced_Renamed)) {
    GCST4023_Reduced_Renamed$eaf.outcome <- rep(NA, nrow(GCST4023_Reduced_Renamed))
  }
  
  # Harmonise Data between exposure and outcome datasets
  harmonised <- harmonise_data(exposure_dat = eQTL_WB, outcome_dat = GCST4023_Reduced_Renamed, action = 2)
  
  # Save harmonised data
  if (!dir.exists('eQTLgen/MAF/Data_MR/LCL/Harmonised_Data/')) {
    dir.create('eQTLgen/MAF/Data_MR/LCL/Harmonised_Data/', recursive = TRUE)
  }
  saveRDS(harmonised, sprintf('eQTLgen/MAF/Data_MR/LCL/Harmonised_Data/Har_Chr%s.rds', i))
  
  # Perform Mendelian Randomization
  if (!dir.exists('eQTLgen/MAF/Data_MR/LCL/MR_Result/')) {
    dir.create('eQTLgen/MAF/Data_MR/LCL/MR_Result/', recursive = TRUE)
  }
  result <- mr(harmonised)
  saveRDS(result, sprintf('eQTLgen/MAF/Data_MR/LCL/MR_Result/Raw_MR_Chr%s.rds', i))
}

# Initialize with the expected columns (adjust these as needed)
result_signif_total <- data.frame(
  id.exposure = character(),
  pval = numeric(),
  chromosome = character(),
  stringsAsFactors = FALSE
)

# Loop to read and append rows
#seuil <- 0.05 / length(unique(clumped_exposure01$id.exposure))  
seuil <- 1
for (i in 1:22) {
  file_path <- sprintf('eQTLgen/MAF/Data_MR/LCL/MR_Result/Raw_MR_Chr%s.rds', i)
  
  if (!file.exists(file_path)) next
  
  result_temp <- readRDS(file_path)
  
  # Make sure result_temp has the expected column; adjust filter as necessary.
  if (nrow(result_temp %>% filter(result_temp$pval < seuil)) > 0) {
    result_signif <- result_temp %>% filter(result_temp$pval < seuil)
    result_signif$chromosome <- sprintf("%s", i)
    result_signif_total <- rbind(result_signif_total, result_signif)
  }
}

# Check structure after reading in the data
print(str(result_signif_total))

# Now, before looping for gene names, ensure result_signif_total has the required column:
if (!"id.exposure" %in% colnames(result_signif_total)) {
  stop("result_signif_total does not contain the column 'id.exposure'")
}

# Process each row to add gene names
rownames(result_signif_total) <- NULL  # Reset row names
result_signif_total_WithNames <- data.frame()

gtf_reduced <- readRDS('GeneCode_ConversionTable.rds')

for (j in 1:nrow(result_signif_total)) {
  row_j <- result_signif_total[j, , drop = FALSE]
  
  # Check that the row is not empty and that id.exposure exists
  if (ncol(row_j) == 0) {
    message("Row ", j, " has no columns. Skipping.")
    next
  }
  
  id_exp <- row_j$id.exposure
  if (length(id_exp) == 0) {
    id_exp <- ""
  }
  id_exp <- as.character(id_exp)
  
  if (nchar(id_exp) > 0) {
    parts <- unlist(strsplit(id_exp, split = "\\."))
    id_part <- if (length(parts) >= 1) parts[1] else id_exp
  } else {
    id_part <- id_exp
  }
  
  entry <- gtf_reduced[gtf_reduced$id == id_part, ]
  if (nrow(entry) == 0 || length(entry$name) == 0) {
    message("Not found for id: ", id_exp)
    gene_name <- id_exp
  } else {
    gene_name <- entry$name[1]
  }
  
  temp_name <- cbind(row_j, geneName = gene_name)
  result_signif_total_WithNames <- rbind(result_signif_total_WithNames, temp_name)
  message("found")
}

# Save or further process result_signif_total_WithNames as needed

names(result_signif_total_WithNames)[names(result_signif_total_WithNames) == "name"] <- "geneName"
result_signif_total <- result_signif_total_WithNames
saveRDS(result_signif_total, 'eQTLgen/MAF/Data_MR/LCL/MR_Result/allResultWithNames.rds')

# Checking Replication ---------------------------------------------------------
#seuil_eQTL <- 0.05 / 50
#MR4 <- readRDS('eQTLgen/MAF/Data_MR/LCLas/MR_Result/allResultSignifWithNames.rds')
#MR4_filtered <- MR4[MR4$pval < seuil_eQTL, ]
write.table(unique(result_signif_total$geneName), 'eQTLgen/LCLGenes_nonsignif.txt', 
          quote = FALSE, row.names = FALSE, col.names = FALSE)
#saveRDS(result_signif_total, 'eQTLgen/MAF/Data_MR/LCLas/MR_Result/allResultSignifWithNames.rds')

# Initialize with the expected columns (adjust these as needed)
result_signif_total <- data.frame(
  id.exposure = character(),
  pval = numeric(),
  chromosome = character(),
  stringsAsFactors = FALSE
)

seuil <- 0.05 / length(unique(clumped_exposure01$id.exposure))  
#seuil <- 1
for (i in 1:22) {
  file_path <- sprintf('eQTLgen/MAF/Data_MR/LCL/MR_Result/Raw_MR_Chr%s.rds', i)
  
  if (!file.exists(file_path)) next
  
  result_temp <- readRDS(file_path)
  
  # Make sure result_temp has the expected column; adjust filter as necessary.
  if (nrow(result_temp %>% filter(result_temp$pval < seuil)) > 0) {
    result_signif <- result_temp %>% filter(result_temp$pval < seuil)
    result_signif$chromosome <- sprintf("%s", i)
    result_signif_total <- rbind(result_signif_total, result_signif)
  }
}

# Check structure after reading in the data
print(str(result_signif_total))

# Now, before looping for gene names, ensure result_signif_total has the required column:
if (!"id.exposure" %in% colnames(result_signif_total)) {
  stop("result_signif_total does not contain the column 'id.exposure'")
}

# Process each row to add gene names
rownames(result_signif_total) <- NULL  # Reset row names
result_signif_total_WithNames <- data.frame()

gtf_reduced <- readRDS('GeneCode_ConversionTable.rds')

for (j in 1:nrow(result_signif_total)) {
  row_j <- result_signif_total[j, , drop = FALSE]
  
  # Check that the row is not empty and that id.exposure exists
  if (ncol(row_j) == 0) {
    message("Row ", j, " has no columns. Skipping.")
    next
  }
  
  id_exp <- row_j$id.exposure
  if (length(id_exp) == 0) {
    id_exp <- ""
  }
  id_exp <- as.character(id_exp)
  
  if (nchar(id_exp) > 0) {
    parts <- unlist(strsplit(id_exp, split = "\\."))
    id_part <- if (length(parts) >= 1) parts[1] else id_exp
  } else {
    id_part <- id_exp
  }
  
  entry <- gtf_reduced[gtf_reduced$id == id_part, ]
  if (nrow(entry) == 0 || length(entry$name) == 0) {
    message("Not found for id: ", id_exp)
    gene_name <- id_exp
  } else {
    gene_name <- entry$name[1]
  }
  
  temp_name <- cbind(row_j, geneName = gene_name)
  result_signif_total_WithNames <- rbind(result_signif_total_WithNames, temp_name)
  message("found")
}

# Save or further process result_signif_total_WithNames as needed

names(result_signif_total_WithNames)[names(result_signif_total_WithNames) == "name"] <- "geneName"
result_signif_total <- result_signif_total_WithNames
saveRDS(result_signif_total, 'eQTLgen/MAF/Data_MR/LCL/MR_Result/allResultSignifWithNames.rds')

# Checking Replication ---------------------------------------------------------
#seuil_eQTL <- 0.05 / 50
#MR4 <- readRDS('eQTLgen/MAF/Data_MR/LCLas/MR_Result/allResultSignifWithNames.rds')
#MR4_filtered <- MR4[MR4$pval < seuil_eQTL, ]
write.table(unique(result_signif_total$geneName), 'eQTLgen/LCLGenes_signif.txt', 
          quote = FALSE, row.names = FALSE, col.names = FALSE)
#saveRDS(result_signif_total, 'eQTLgen/MAF/Data_MR/LCLas/MR_Result/allResultSignifWithNames.rds')
