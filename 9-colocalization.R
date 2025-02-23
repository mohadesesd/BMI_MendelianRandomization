# Header ----------------------------------------------------------------------

# AUTHOR: John Eric Hamilton
# PURPOSE: Performing colocalization study (inspired by code from Basile 02 BJ Colocalization),
#          looping through all significant genes
# NOTES: 
#   *Inputs:   - MR Results filtered and non-filtered (to choose one gene at a time)
#              - Exposition file 
#              - GWAS (Chiou Only)
#
#   *Outputs:  - Colocalization test statistics  
#              - SNP plot with p-value

# Loading Library and Data -----------------------------------------------------
library(arrow)
library(coloc)
library(susieR)
library(gridExtra)
library(tidyverse)
library(LDlinkR) #To get the LD
library('reshape') #Formating data
library(TwoSampleMR)
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(reshape2)
library(LDlinkR)
library(ieugwasr)
library(httr)
library(purrr)
library(jsonlite)
# Importing the harmonized data for chromosome 1 to 22 and merging all
Har_total <- data.frame()
for (i in 1:22){
  result_temp <- readRDS(sprintf('./Pancreas_MR/Har/Har_Chr%s.rds', i))
  result_temp$chromosome <- sprintf("%s", i)
  Har_total <- rbind(Har_total, result_temp)
}
# Initialize result table
tab_res <- data.frame()

# Significant results that were not in MHC region
signif_total <- readRDS('./Pancreas_MR/significant_allResultSignifWithNames.rds')
merged_df <- inner_join(Har_total, signif_total, by = c("id.exposure" = "id.exposure"))
exposure <- read.delim("/lustre03/project/6050805/SHARED_FULL/SHARED/GTEx/sQTL_GTEx_v8/GTEx_Analysis_v8_sQTL_independent/Pancreas.v8.independent_sqtls.txt.gz")
# Loop -------------------------------------------------------------------------
# Taking the gene TH to start (eventually will loop through all of them)
for (i in 1:length(merged_df$geneName)){
  # Gene to coloc analyze information
  ensembleID <- sub(".*:(ENSG[0-9]+)\\..*", "\\1", merged_df$id.exposure)[i]
  print(ensembleID)
  prot <- merged_df$geneName[i]
  print(prot)
  entry <- merged_df[grep(ensembleID, merged_df$phenotype_id),]
  chrI <- colsplit(entry$variant_id.y, "_", c("chromosome", "base_pair_location", "other_allele", "effect_allele", "name"))[1][i]
  position <-colsplit(entry$variant_id.y, "_", c("chromosome", "base_pair_location", "other_allele", "effect_allele", "name"))[2][i]
  SNPI <- paste(chrI, ":", position[1,1], sep="")
  print(SNPI)
  # 1000000 Mb windows
  pos_min <- position[1,1] - 500000
  pos_max <- position[1,1] + 500000
  
  ### EXPO FILE ###

  pattern <- sprintf("^%s_", chrI)
  #print(pattern)
  # Pre-processing exposure
  ## Getting a single chromosome
  sQTL_Artery <- exposure %>% 
    mutate(variant_id = as.character(variant_id)) %>% 
    filter(str_detect(variant_id, pattern))
  ## Getting the columns of interest
  sQTL_Artery_Reduced <- select(sQTL_Artery, phenotype_id, variant_id, maf, pval_nominal, slope, slope_se, ref_factor, tss_distance, ma_count, ma_samples)
  sQTL_Artery_Reduced$chromosome <- colsplit(sQTL_Artery_Reduced$variant_id, "_", 
    c("chromosome", "base_pair_location", "other_allele", "effect_allele", "name"))[[1]]
  sQTL_Artery_Reduced$base_pair_location <- colsplit(sQTL_Artery_Reduced$variant_id, "_", 
    c("chromosome", "base_pair_location", "other_allele", "effect_allele", "name"))[[2]]
  sQTL_Artery_Reduced$other_allele <- colsplit(sQTL_Artery_Reduced$variant_id, "_", 
    c("chromosome", "base_pair_location", "other_allele", "effect_allele", "name"))[[3]]
  sQTL_Artery_Reduced$effect_allele <- colsplit(sQTL_Artery_Reduced$variant_id, "_", 
    c("chromosome", "base_pair_location", "other_allele", "effect_allele", "name"))[[4]]
  sQTL_Artery_Reduced$SNP <- paste(sQTL_Artery_Reduced$chromosome, ":", sQTL_Artery_Reduced$base_pair_location, sep="")

  # Changes the MAF based on the value of ref_factor (either -1 or 1)
  sQTL_Artery_Reduced[sQTL_Artery_Reduced$ref_factor == -1,]$maf <- 1- sQTL_Artery_Reduced[sQTL_Artery_Reduced$ref_factor == -1,]$maf  
  # cut the Exposure for the given base pair range
  expo_cut <- sQTL_Artery_Reduced[sQTL_Artery_Reduced$base_pair_location <= pos_max & sQTL_Artery_Reduced$base_pair_location >= pos_min, ]
  expo_cut = expo_cut[expo_cut$other_allele%in%c('A','T','C','G') & expo_cut$effect_allele%in%c('A','T','C','G'),]
  expo_cut <- cbind(expo_cut[,c(15,1,11,12,13,14,8,10,9,3,4,5,6)])
   ## Remove Palindromic Variants
  #expo_cut <- expo_cut[
   # !((expo_cut$effect_allele.exposure == "A" & expo_cut$other_allele.exposure == "T") |
    #  (expo_cut$effect_allele.exposure == "T" & expo_cut$other_allele.exposure == "A") |
     # (expo_cut$effect_allele.exposure == "C" & expo_cut$other_allele.exposure == "G") |
     # (expo_cut$effect_allele.exposure == "G" & expo_cut$other_allele.exposure == "C")), ]

  ### Renaming columns for harmonize_data() function
  colnames(expo_cut) <- c("SNP", "phenotype_id", "chromosome", "base_pair_location", "other_allele.exposure", "effect_allele.exposure", "tss_distance", "ma_samples", "ma_count", "eaf.exposure", "pval.exposure", "beta.exposure", "se.exposure")
  print(expo_cut$se.exposure)
  ### OUTCOME FILE ###
  outc_interest <- readRDS(sprintf('./Pancreas_MR/GCST2409/GCST2409_%s.rds', chrI))
  print(head(outc_interest))
  proxied <- readRDS('./Pancreas_MR/outcomeToMerge.rds')
  outc_interest <- rbind(outc_interest, proxied %>% filter(sprintf("%s", chrI) == proxied$chromosome))
  print(head(outc_interest))
  
  
  # Remove Indel
  outc_interest=outc_interest[outc_interest$other_allele%in%c('A','T','C','G') & outc_interest$effect_allele%in%c('A','T','C','G'),]
  
  # Allele mineur rare (to remove some of the duplicates in the outcome)
  #outc_interest <- outc_interest[outc_interest$effect_allele_frequency > 0.01,]
  
  # merge the 2 GWAS using harmonise_data() function 
  #expo_cut_Renamed <- expo_cut[expo_cut$gene_id == ensembleID,]
  expo_cut$id.exposure <- expo_cut$phenotype_id
  expo_cut$exposure <- expo_cut$phenotype_id
    if (!"eaf.outcome" %in% names(outc_interest)) {
  outc_interest$eaf.outcome <- rep(NA, nrow(outc_interest))
  }
  
  harm_data <- harmonise_data(expo_cut, outc_interest, action=2)
  saveRDS(harm_data, sprintf('./Pancreas_MR/Genes/Coloc_Harm_%s_result.rds', prot))
  
  
  # Remove duplicated SNPs that still remain
  test <- duplicated(harm_data$SNP)
  SNP_toRemove <- harm_data[test,"SNP"]
  harm_data <- harm_data[!(harm_data$SNP %in% SNP_toRemove),]
  
  # Rearrange data
  tmp <- harm_data[, c(1, 29, 16, 17, 6, 28, 27, 25, 8, 7, 18, 15)]
  
  colnames(tmp) <- c("SNP", "gene_id", "chr", "position", "expo_eff", "expo_se", "expo_pv", "expo_N", "expo_maf",  
                     "outc_eff", "outc_se", "outc_pv")
  # varbeta = square(se)
  tmp$expo_varbeta <- (tmp$expo_se)^2
  tmp$outc_varbeta <- (tmp$outc_se)^2
  
  # COLOC ANALYSIS
  # 2 dataset for analyse Prot - ANM
  D1_expo <- list() 
  D2_outc <- list()
  
  # expo
  D1_expo$type <- "quant"
  D1_expo$snp <- tmp$SNP
  D1_expo$position <- tmp$position
  D1_expo$beta <- tmp$expo_eff
  D1_expo$varbeta <- tmp$expo_varbeta
  D1_expo$varbeta[is.na(D1_expo$varbeta)] <- 1
  D1_expo$sdY <- 1 # (prot level are scale)
  
  
  # outc
  D2_outc$type <- "cc"
  D2_outc$snp <- tmp$SNP
  D2_outc$position <- tmp$position
  D2_outc$beta <- tmp$outc_eff
  D2_outc$varbeta <- tmp$outc_varbeta
  D2_outc$sdY <- 1 # (FROM GWAS STUDY : 1.3 for AAM and 4 for ANM)
  
  
  saveRDS(D1_expo, sprintf('./Pancreas_MR/Genes/Coloc_Expo_%s_result.rds', prot))
  saveRDS(D2_outc, sprintf('./Pancreas_MR/Genes/Coloc_Outc_%s_result.rds', prot))
  
  # Run coloc analysis 
  res <- coloc.abf(dataset1 = D1_expo, dataset2 = D2_outc)
  
  # save the result
  tab_res <- rbind(tab_res, c(paste(prot), as.vector(res$summary)))
  
}

# Format and Save Results
colnames(tab_res) <- c("Gene Name", "nSNP", "H0", "H1", "H2", "H3", "H4")
tab_res$H0 <- as.numeric(tab_res$H0)
tab_res$H1 <- as.numeric(tab_res$H1)
tab_res$H2 <- as.numeric(tab_res$H2)
tab_res$H3 <- as.numeric(tab_res$H3)
tab_res$H4 <- as.numeric(tab_res$H4)

saveRDS(tab_res, './Pancreas_MR/Coloc_results/Coloc_result.rds')

# Filter for genes that passed colocalization 
final_res <- tab_res[tab_res$H4 > 0.8,] 
final_moreAnalysis <- tab_res[!(tab_res$H4 > 0.8),]

saveRDS(final_res,'./Pancreas_MR/Coloc_results/passingColocGenes.rds')
saveRDS(final_moreAnalysis,'./Pancreas_MR/Coloc_results/ColocGenesToProxy.rds')


write.table(final_res$`Gene Name`, './Pancreas_MR/Coloc_results/PassingGeneNames.txt', quote = FALSE, row.names = FALSE, col.names = FALSE)

final_res <- readRDS('./Pancreas_MR/Coloc_results/passingColocGenes.rds')
final_moreAnalysis <- readRDS('./Pancreas_MR/Coloc_results/ColocGenesToProxy.rds')
final_all <- rbind(final_res, final_moreAnalysis)

#Filter for genes that might need proxy (High H3 and very low H4)
final_moreAnalysis_Proxy <- final_moreAnalysis[final_moreAnalysis$H3 > 0.8,]
