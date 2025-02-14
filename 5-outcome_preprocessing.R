# Header ---------------------------------------------------------------------

# AUTHOR: John Eric Hamilton
# PURPOSE: Preprocessing outcome data
# NOTES: 
#   *Inputs:   - Outcome data 
#
#   *Outputs:  - Preprocessed_Data

# Load libraries and functions -------------------------------------------------
library(TwoSampleMR)
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(reshape2)

freadf <- function(df){return(as.data.frame(fread(df)))}

# Load datasets ----------------------------------------------------------------
df_GCST90002409 <- freadf("GCST90002409_buildGRCh38.tsv")
print(head(df_GCST90002409))

# Run MR -----------------------------------------------------------------------
for (i in 1:22){
  # Pre-processing outcome
  ## Getting a single chromosome
  GCST2409 <- df_GCST90002409 %>% filter(sprintf("chr%s", i) == df_GCST90002409$chromosome)
  
  ### To upper case
  GCST2409$effect_allele <- toupper(GCST2409$effect_allele)
  GCST2409$other_allele  <- toupper(GCST2409$other_allele)
  print(head(GCST2409))
  ### Making SNP column
  GCST2409$SNP <- paste(GCST2409$chromosome, ":", GCST2409$base_pair_location, sep="")
  ## Merging Proxied Data
  #proxied <- readRDS('outcomeToMerge.rds')
  #proxied <- select(proxied, -SNP_Original)
  #GCST2409 <- rbind(GCST2409, proxied %>% filter(sprintf("%s", i) == proxied$chromosome))
  
  ## Getting the columns of interest
  #GCST2409_Reduced <- select(GCST2409, -sample_size)
  
  # Remove Indel
  GCST2409=GCST2409[GCST2409$other_allele%in%c('A','T','C','G') & GCST2409$effect_allele%in%c('A','T','C','G'),]

  ### Reduced and ordered table

  GCST2409_Renamed <- cbind(GCST2409[,c(1,10,5,4,8,6,7,9)])
  ### Renaming columns for harmonize_data() function
  colnames(GCST2409_Renamed) <- c("variant_id", "SNP", "other_allele.outcome",  "effect_allele.outcome", "pval.outcome", "beta.outcome", "se.outcome", 'samplesize.outcome')
  
  # Remove palindromic variants
  GCST2409_Renamed <- GCST2409_Renamed[
    !((GCST2409_Renamed$effect_allele.outcome == "A" & GCST2409_Renamed$other_allele.outcome == "T") |
      (GCST2409_Renamed$effect_allele.outcome == "T" & GCST2409_Renamed$other_allele.outcome == "A") |
      (GCST2409_Renamed$effect_allele.outcome == "C" & GCST2409_Renamed$other_allele.outcome == "G") |
      (GCST2409_Renamed$effect_allele.outcome == "G" & GCST2409_Renamed$other_allele.outcome == "C")), ]
  
  ### Adding outcome column
  GCST2409_Renamed$outcome <- "BMI" 
  GCST2409_Renamed$id.outcome <- "BMI"
  print(head(GCST2409_Renamed))

  saveRDS(GCST2409_Renamed, sprintf('GCST2409_chr%s.rds', i))
}
