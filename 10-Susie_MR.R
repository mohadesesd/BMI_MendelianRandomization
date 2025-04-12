# Header ----------------------------------------------------------------------

# AUTHOR: John Eric Hamilton
# PURPOSE: Performing Susie colocalization study 
#
# NOTES: 
#   *Inputs:   - LD Matrix
#              - Harm Data from Exposition + Outcome of Chiou only
#              - Plink BED file
#
#   *Outputs:  - Colocalization statistics  
#              
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
# Prepare Data to be used for co localization (RSID) ---------------------------  
Har_total <- data.frame()
for (i in 1:22){
  result_temp <- readRDS(sprintf('./Pituitary_MR/Har/Har_Chr%s.rds', i))
  result_temp$chromosome <- sprintf("%s", i)
  Har_total <- rbind(Har_total, result_temp)
}

# Significant results that were not in MHC region (That didnt pass coloc that we want to do Susie Coloc with)
signif_total <- readRDS('./Pituitary_MR/significant_allResultSignifWithNames.rds')
merged_df <- inner_join(Har_total, signif_total, by = c("id.exposure" = "id.exposure"))
final_moreAnalysis <- readRDS('./Pituitary_MR/Coloc_results/ColocGenesToProxy.rds')
merged_df <- merged_df[(merged_df$geneName %in% final_moreAnalysis$`Gene Name`),]
exposure <- read.delim("/lustre03/project/6050805/SHARED_FULL/SHARED/GTEx/sQTL_GTEx_v8/GTEx_Analysis_v8_sQTL_independent/Pituitarycutaneous.v8.independent_sqtls.txt.gz")

# To generate the new tmp files to be used with Plink
for (i in 1:length(merged_df$geneName)){
  
  # Get information of gene of interest
  # Taking the gene ADCY3 to start (eventually will loop through all of them)
  print(i)
  ensembleID <- sub(".*:(ENSG[0-9]+)\\..*", "\\1", merged_df$id.exposure)[i]
  print(ensembleID)
  prot <- merged_df$geneName[i]
  print(prot)
  entry <- merged_df[grep(ensembleID, merged_df$phenotype_id),]
  split_data <- colsplit(entry$variant_id.y, "_", 
                       c("chromosome", "base_pair_location", "other_allele", "effect_allele", "name"))
  chrI <- split_data[1, "chromosome"]
  position <- split_data[1, "base_pair_location"]
  print(position)
  SNPI <- paste(chrI, ":", position, sep="")
  print(SNPI)
  
  # 1000000 Mb windows
  pos_min <- position - 500000
  pos_max <- position + 500000
  
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
  colnames(expo_cut) <- c("SNP", "phenotype_id", "chromosome", "base_pair_location", "other_allele.exposure", "effect_allele.exposure", "tss_distance", "ma_samples", "ma_count", "eaf.exposure", "pval.exposure", "beta.exposure", "se.exposure")

   ## Remove Palindromic Variants
  expo_cut <- expo_cut[
   !((expo_cut$effect_allele.exposure == "A" & expo_cut$other_allele.exposure == "T") |
    (expo_cut$effect_allele.exposure == "T" & expo_cut$other_allele.exposure == "A") |
    (expo_cut$effect_allele.exposure == "C" & expo_cut$other_allele.exposure == "G") |
    (expo_cut$effect_allele.exposure == "G" & expo_cut$other_allele.exposure == "C")), ]
  ### OUTCOME FILE ###
  outc_interest <- readRDS(sprintf('./Pituitary_MR/GCST2409/GCST2409_%s.rds', chrI))
  
  expo_cut$id.exposure <- expo_cut$phenotype_id
  expo_cut$exposure <- expo_cut$phenotype_id
  expo_cut <- expo_cut[sub(".*:(ENSG[0-9]+)\\..*", "\\1", expo_cut$id.exposure) == ensembleID,]
  if (!"eaf.outcome" %in% names(outc_interest)) {
  outc_interest$eaf.outcome <- rep(NA, nrow(outc_interest))
  }
    
  outc_interest_Renamed <- outc_interest
  outc_interest_Renamed$id.outcome <- "BMI"
  outc_interest_Renamed$outcome <- "BMI"
  
  harm_data <- harmonise_data(expo_cut, outc_interest_Renamed)
  
  # Remove duplicated SNPs that still remain
  test <- duplicated(harm_data$SNP)
  SNP_toRemove <- harm_data[test,"SNP"]
  harm_data <- harm_data[!(harm_data$SNP %in% SNP_toRemove),]
  # rearrange data
  tmp <- harm_data[, c(1, 28, 21, 22, 6, 26, 25, 23, 8, 7, 16, 15)]
  
  colnames(tmp) <- c("SNP", "gene_id", "chr", "position", "expo_eff", "expo_se", "expo_pv", "expo_N", "expo_maf",  
                     "outc_eff", "outc_se", "outc_pv")
  
  # varbeta = square(se)
  tmp$expo_varbeta <- (tmp$expo_se)^2
  tmp$outc_varbeta <- (tmp$outc_se)^2
  ### AT THIS POINT TMP FILES ARE GENERATED ###
  
  # Get RS ID using outcome file and add that column in tmp
  outc_RSID <- outc_interest  
  # Identify Similar SNPs 
  outc_RSID_row <- outc_RSID[outc_RSID$SNP%in%tmp$SNP,]
  outc_RSID_row <- outc_RSID_row %>% select("variant_id", "SNP")
  outc_RSID_row <- unique(outc_RSID_row)
  # Merge by SNP to add the RSID column to the data file 
  tmp_new <- merge(harm_data, outc_RSID_row, by.x="SNP", by.y="SNP")
  print(ncol(tmp_new))
  saveRDS(tmp_new, sprintf("./Pituitary_MR/Plink8/tmp/tmp_%s_%s.rds", chrI, prot))
  write.table(tmp_new$variant_id.x, sprintf('./Pituitary_MR/Plink8/test%s_%s.txt', chrI, prot), row.names = FALSE, col.names = FALSE, quote = FALSE) #Use in Plink
}
#plink --bfile "~/projects/def-dman2020/SHARED_FULL/SHARED/1000Genomes/wgs/plink/EUR/Chr2.EUR" --extract "./Pituitary_MR/Plink8/testchr2_ADCY3.txt" --recodeA --out "susieChr2"


merged_df <- merged_df[!(duplicated(merged_df$geneName)),]
susieAllResults <- NULL

for (i in 1:length(merged_df$geneName)){
  print(i)
  ensembleID <- sub(".*:(ENSG[0-9]+)\\..*", "\\1", merged_df$id.exposure)[i]
  print(ensembleID)
  prot <- merged_df$geneName[i]
  print(prot)
  entry <- merged_df[grep(ensembleID, merged_df$phenotype_id),]
  split_data <- colsplit(entry$variant_id.y, "_", 
                       c("chromosome", "base_pair_location", "other_allele", "effect_allele", "name"))
  chrI <- split_data[1, "chromosome"]
  position <- split_data[1, "base_pair_location"]
  print(position)
  SNPI <- paste(chrI, ":", position, sep="")
  print(SNPI)
  
  # And after estimate LD with R
  plinkLD_raw <- function(raw_files = sprintf("./Pituitary_MR/Plink8/susie%s_%s.raw", chrI, prot)) {
    #' fn : imput mean for NA data
    #'
    #' @param x a matrix
    #'
    #' @export
    fn <- function(x){
      m. <- mean(x, na.rm = TRUE)
      x[is.na(x)] <- m.
      return(x)
    }
    t1 <- read.table(raw_files, header = T)
    ld1 <- as.matrix(t1[, c(7:ncol(t1))])
    # number of NA
    nNA <- sum(is.na(ld1))
    nNA <- nNA / length(as.vector(ld1)) * 100
    # remplace NA by colmeans
    ld1 <- fn(ld1)
    # add intercept (just a trick to avaid null ecart type)
    ld1 <- rbind(rep(0.5, ncol(ld1)),
                 ld1)
    # R value (correlation MATRIX)
    ld1 <- cor(ld1)
    # issues with the rsid (need remove _A, _T, ...)
    rs <- data.frame(rs = colnames(ld1))
    rs <- tidyr::separate(rs, rs, c("rs", "other"), "_")
    colnames(ld1) <- rs$rs
    row.names(ld1) <- rs$rs
    return(list(ld_r = ld1,
                out_plink = t1,
                nNA = paste0(round(nNA, 4), " % of NAs")))
  }
  
  # Define Function to remove NA values
  fn <- function(x){
    
    m. <- mean(x, na.rm = TRUE)
    x[is.na(x)] <- m.
    return(x)
  }
  
  # Loading the file from plink into RStudio
  tmps <- plinkLD_raw(sprintf("./Pituitary_MR/Plink8/susie%s_%s.raw", chrI, prot))
  
  
  # Making the LD matrix ---------------------------------------------------------
  ld_r <- tmps$ld_r
  
  # For outc_plink
  outc_plink <- tmps$out_plink
  outc_plink <- data.frame(rsid = colnames(outc_plink)[7:ncol(outc_plink)])
  outc_plink <- separate(outc_plink, rsid, c("rsid", "eff_all"), "_")
  colnames(outc_plink) <- c("rsid", "ld_eff_all")
  outc_plink <- outc_plink[outc_plink$rsid %in% colnames(tmps$ld_r), ]
  
  tmp_new <- readRDS(sprintf("./Pituitary_MR/Plink8/tmp/tmp_%s_%s.rds", chrI, prot))
  # Merge outc_plink with tmp file that has the RSIDs
  tmp_new_susie <- merge.data.frame(outc_plink, tmp_new, by.x = 1, by.y=32)
  
  # remove mismatch between eff allele expo, ouctome and ld matrix (metabol = effect_allele)
  dat_susie <- tmp_new_susie[(tmp_new_susie$ld_eff_all == tmp_new_susie$effect_allele.exposure) & (tmp_new_susie$ld_eff_all ==   tmp_new_susie$effect_allele.outcome), ]
  
  # ld with the right set of rsid
  ld_r <- ld_r[dat_susie$rsid, dat_susie$rsid] #LD MATRIX
  
  saveRDS(dat_susie, sprintf('./Pituitary_MR/Plink8/susieResults/dat_susie%s_%s.rds', chrI, prot))
  saveRDS(ld_r, sprintf('./Pituitary_MR/Plink8/susieResults/LD_MatrixChr%s_%s.rds', chrI, prot))
  
  # calculate Z-score for exposure data and outcome data
  dat_susie$z_expo <- dat_susie$beta.exposure / dat_susie$se.exposure
  dat_susie$z_outc <- dat_susie$beta.outcome / dat_susie$se.outcome
  print(head(dat_susie))
  colnames(dat_susie) <- c("rsid", "ld_eff_all", "snp", "effect_allele", "other_allele", "effect_allele.outcome", 
  "other_allele.outcome", "beta.exposure", "beta.outcome", "eaf.exposure", "eaf.outcome", "remove", "palindromic", 
  "ambiguous", "id.outcome","variant_id.x", "pval.outcome", "se.outcome",
  "sample_size", "outcome", "phenotype_id", "chromosome", "base_pair_location", "tss_distance", "ma_samples", 
  "ma_count", "pval.exposure",  "se.exposure",
   "id.exposure", "exposure", "action", "SNP_index", "mr_keep", "z_expo", "z_outc")
  
  # Check for NA
  ld_r_na <- fn(ld_r)
  
  # calculate lambda
  # concordance between the LD matrix and the z-scores (should be near 0)
  lambda_expo <- estimate_s_rss(dat_susie$z_expo, (ld_r_na)^2, n = mean(dat_susie$ma_samples))
  lambda_outc <- estimate_s_rss(dat_susie$z_outc, (ld_r_na)^2, n = mean(dat_susie$sample_size))
  
  
  
  # use susis_rss same as runsusie (runsusie add some names)
  fit_expo <- susie_rss(dat_susie$z_expo, (ld_r_na)^2, n = mean(dat_susie$ma_samples))
  fit_outc <- susie_rss(dat_susie$z_outc, (ld_r_na)^2, n = mean(dat_susie$sample_size))
  
  
  # Coloc susie ------------------------------------------------------------------
  res <- coloc.susie(fit_expo, fit_outc) 
  
  #View(res$summary)
  
  
  saveRDS(res, sprintf("./Pituitary_MR/Plink8/susieResults/coloc_susie_Chr%s_%s.rds", chrI, prot))
  saveRDS(res$summary , sprintf("./Pituitary_MR/Plink8/susieResults/susieResults_%s_%s.rds", chrI, prot))
  
  susieSignifResults <- res$summary[res$summary$PP.H4.abf > 0.8, ]
  susieAllResults <- rbind(susieAllResults, susieSignifResults)
  print("Next one...")
}

# Summary Table
susieAllResults_ <- data.frame()
for (i in 1:length(signif_total$geneName)){
  ensembleID <- signif_total$id.exposure[i]
  prot <- signif_total$geneName[i]
  chrI <- signif_total$chromosome[i]
  entry <- df[grep(ensembleID, df$gene_id),]
  position <-colsplit(entry$variant_id, "_", c("chromosome", "base_pair_location", "other_allele", "effect_allele", "name"))[2]
  SNPI <- paste0("chr", chrI, ":", position[1,1])
  
  
  res <- readRDS(sprintf("./Pituitary_MR/Plink8/susieResults/coloc_susie_Chr%s_%s.rds", chrI, prot))
  
  if (length(res) != 1){
    susieSignifResults <- res$summary[res$summary$PP.H4.abf > 0.8, ]
    susieSignifResults$Gene <- paste(prot)
    susieAllResults_ <- rbind(susieAllResults_, susieSignifResults)
  }
}

saveRDS(susieAllResults_, "./Pituitary_MR/Plink8/susieResults/allSignificantSusiesChrome2_6.rds")
