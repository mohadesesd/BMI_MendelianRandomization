# Header ---------------------------------------------------------------------

# AUTHOR: John Eric Hamilton
# PURPOSE: Interpretting the results of MR
# NOTES: 
#   *Inputs:   - MR Results 
#              - Seuil (Bonferroni calculated by dividing 0.05 by the number of entries in the exposure)
#
#   *Outputs:  - All significant resulting SNPs from the MR (called result_signif_total)
#              - HLA/MHC region filtered significant SNPs (called filtered_MHC)
#              - non-HLA region filtered significant SNPs (called filtered) 
#   *Resulting files contain merged chromosome results

library(TwoSampleMR)
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(reshape2)
library(LDlinkR)
# Generating Results ----------------------------------------------------------

# Recall the Bonferroni threshold (sQTL Artery)
seuil <- 0.05

# Importing the results for chromosome 1 to 22 and merging all significant genes in one table
result_signif_total <- data.frame()
result_all <- data.frame()
for (i in 1:22){
  result_temp <- readRDS(sprintf('./Pancreas_MR/Raw_MR/Raw_MR_Chr%s.rds', i))
  if (dim(result_temp %>% filter(result_temp$pval < seuil))[1] > 0){
    result_signif <- result_temp %>% filter(result_temp$pval < seuil)
    result_signif$chromosome <- sprintf("%s", i)
    
    result_signif_total <- rbind(result_signif_total, result_signif)
  }
  
  result_temp$chromosome <- sprintf("%s", i)
  result_all <- rbind(result_all, result_temp)
}

saveRDS(result_all,'./Pancreas_MR/allResult_allNonSignifAndSignif.rds')
print(nrow(result_all))
print(nrow(result_signif_total))

## Filtering genes that are in the MHC region (chromosome 6)
filtered <- data.frame()
filtered_MHC <- data.frame()

for (j in 1:nrow(result_signif_total)) {
  #for (j in 1:nrow(result_signif_total)) {
  current_exposure <- result_signif_total[["id.exposure"]][j]
  
  # If current_exposure is NA or empty, skip this iteration
  if (length(current_exposure) == 0 || is.na(current_exposure) || current_exposure == "") {
    warning(paste("Skipping row", j, "because current_exposure is empty"))
    next
  }
  
  temp <- result_all %>% filter(result_all$id.exposure == current_exposure)
  chromosome <- colsplit(temp$id.exposure, ":", c("chromosome", "start", "end", "name"))[[1]] 
  location <- colsplit(temp$id.exposure, ":", c("chromosome", "start", "end", "name"))[[3]]

  if (length(location) > 1){
    location <- location[1]
  }
  
  if (length(chromosome) > 1){
    chromosome <- chromosome[1]
  }
  
  if (chromosome == "chr6"){
    if (between(location, 28510120, 33480577)){
      filtered_MHC <- rbind(filtered_MHC, temp)
    } else {
      # Could be chromosome 6 but not MHC region
      filtered <- rbind(filtered, temp)
    }
  } else {
    filtered <- rbind(filtered, temp)
  }
}


saveRDS(result_signif_total, './Pancreas_MR/allResultSignif.rds')
saveRDS(filtered, './Pancreas_MR/allResults_filtered.rds')
saveRDS(filtered_MHC, './Pancreas_MR/allResults_filtered_MHC_region.rds')
