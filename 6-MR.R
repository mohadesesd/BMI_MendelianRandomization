
# Header ---------------------------------------------------------------------

# AUTHOR: John Eric Hamilton
# PURPOSE: Running the Mendelian Randomisation
# NOTES: 
#   *Inputs:   - Exposure and Outcome data 
#
#   *Outputs:  - harmonized data
#              - MR results

# Load libraries and functions -------------------------------------------------
library(TwoSampleMR)
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(reshape2)
library(LDlinkR)

token <- "23a23730fd0b"

freadf <- function(df){return(as.data.frame(fread(df)))}

# Load datasets ----------------------------------------------------------------
df_GCST90002409 <- freadf("GCST90002409_buildGRCh38.tsv")

#print(head(df_GCST90002409))
df <- read.delim("/lustre03/project/6050805/SHARED_FULL/SHARED/GTEx/sQTL_GTEx_v8/GTEx_Analysis_v8_sQTL_independent/Artery_Coronary.v8.independent_sqtls.txt.gz")


exposure <- df

seuil <- 0.05/length(unique(exposure$phenotype_id))

# Run MR -----------------------------------------------------------------------
for (i in 1:22){
  
  sQTL_Artery <- readRDS(sprintf('sQTL_Artery/sQTL_Artery_chr%s.rds', i))
  GCST2409 <- readRDS(sprintf('GCST2409/GCST2409_chr%s.rds', i))
  # Dummy eaf for outcome
  if (!"eaf.outcome" %in% names(GCST2409)) {
  GCST2409$eaf.outcome <- rep(0.5, nrow(GCST2409))
  }
  print(head(sQTL_Artery))
  # Harmonise Data
  harmonised <- harmonise_data(exposure_dat = sQTL_Artery, outcome_dat = GCST2409, action = 2)
  saveRDS(harmonised, sprintf('Har_Chr%s.rds', i))
  
  # Mendelian Randomization
  result <- mr(harmonised)
  saveRDS(result, sprintf('Raw_MR_Chr%s.rds', i))
  
  result_signif <- result %>% filter(result$pval < seuil) 
  saveRDS(result_signif, sprintf('MR_Chr%s.rds', i))
}
