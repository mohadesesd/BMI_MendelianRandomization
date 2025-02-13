
# Header ---------------------------------------------------------------------

# AUTHOR: John Eric Hamilton
# PURPOSE: Preprocessing the exposure
# NOTES: 
#   *Inputs:   - Exposure data
#
#   *Outputs:  - Preprocessed data

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

#Filtering exposure for the Fstat with a threshold of 10
#Finding R2 for all values
exposure$R2 <- (get_r_from_bsen(exposure$slope,exposure$slope_se,328))^2 # CHANGE FOR EACH TISSUE
#Finding the Fstat
exposure$Fstat <- (exposure$R2/1)/((1-exposure$R2)/(328-1-1))

#Filtering
exposure <- exposure[exposure$Fstat > 10,]
seuil <- 0.05/length(unique(exposure$phenotype_id))

# Run MR -----------------------------------------------------------------------
for (i in 1:22){

  pattern <- sprintf("^chr%s_", i)
  # Pre-processing exposure
  ## Getting a single chromosome
  sQTL_Artery <- exposure %>% 
    mutate(variant_id = as.character(variant_id)) %>% 
    filter(str_detect(variant_id, pattern))
  ## Getting the columns of interest
  sQTL_Artery_Reduced <- select(sQTL_Artery, phenotype_id, variant_id, maf, pval_nominal, slope, slope_se, ref_factor)
  # Changes the MAF based on the value of ref_factor (either -1 or 1)
  sQTL_Artery_Reduced[sQTL_Artery_Reduced$ref_factor == -1,]$maf <- 1- sQTL_Artery_Reduced[sQTL_Artery_Reduced$ref_factor == -1,]$maf


  ## Renaming and separating columns for harmonize_data() function 
  ### variant_id column separation and merging 
  v_id <- colsplit(sQTL_Artery_Reduced$variant_id, "_", c("chromosome", "base_pair_location", "other_allele", "effect_allele", "name"))
  v_id$SNP <- paste(v_id$chromosome, v_id$base_pair_location, sep=":")
  
  test_merge <- cbind(sQTL_Artery_Reduced[1], sQTL_Artery_Reduced[2], v_id[,c(6,3,4)], sQTL_Artery_Reduced[,2:5])
  
  ## Remove Indel
  test_merge = test_merge[test_merge$other_allele%in%c('A','T','C','G') & test_merge$effect_allele%in%c('A','T','C','G'),]
  
  ### Renaming columns for harmonize_data() function
  colnames(test_merge) <- c("phenotype_id", "variant_id", "SNP", "other_allele.exposure", "effect_allele.exposure","eaf.exposure", "pval.exposure", "beta.exposure", "se.exposure")
  sQTL_Artery_Reduced_Renamed <- test_merge

  ## Remove Palindromic Variants
  sQTL_Artery_Reduced_Renamed <- sQTL_Artery_Reduced_Renamed[
    !((sQTL_Artery_Reduced_Renamed$effect_allele.exposure == "A" & sQTL_Artery_Reduced_Renamed$other_allele.exposure == "T") |
      (sQTL_Artery_Reduced_Renamed$effect_allele.exposure == "T" & sQTL_Artery_Reduced_Renamed$other_allele.exposure == "A") |
      (sQTL_Artery_Reduced_Renamed$effect_allele.exposure == "C" & sQTL_Artery_Reduced_Renamed$other_allele.exposure == "G") |
      (sQTL_Artery_Reduced_Renamed$effect_allele.exposure == "G" & sQTL_Artery_Reduced_Renamed$other_allele.exposure == "C")), ]
  
  ## Add penotype_id as exposure_id
  sQTL_Artery_Reduced_Renamed$id.exposure <- sQTL_Artery_Reduced_Renamed$phenotype_id
  sQTL_Artery_Reduced_Renamed$exposure <- sQTL_Artery_Reduced_Renamed$phenotype_id
  print(head(sQTL_Artery_Reduced_Renamed))

  saveRDS(sQTL_Artery_Reduced_Renamed, sprintf('sQTL_Artery_chr%s.rds', i))
}