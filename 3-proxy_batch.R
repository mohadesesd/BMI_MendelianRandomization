# AUTHOR: John Eric Hamilton
# PURPOSE: Proxy Search for SNPs that are in the exposure and are missing the outcome
# NOTES: 
#   *Inputs:   - exposure
#              - outcome
#
#   *Outputs:  - Proxy Data 
# Load library -----------------------------------------------------------------
library("LDlinkR")
library(TwoSampleMR)
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(reshape2)
# Loading Data and looking for missing SNPs in outcome ---------------------------------
### GWAS  'df_GCST90002409' (outcome) vs sQTL 'df' (exposure)
freadf <- function(df){return(as.data.frame(fread(df)))}

# Load datasets ----------------------------------------------------------------
df_GCST90002409 <- freadf("GCST90002409_buildGRCh38.tsv")
df_GCST90002409$SNP <- paste("chr", df_GCST90002409$chromosome, ":", df_GCST90002409$base_pair_location, sep="")

sQTL <- data.frame()
for (i in 1:22){
  sQTL <- rbind(sQTL, readRDS(sprintf('sQTL_Artery/sQTL_Artery_chr%s.rds', i)))
}
proxy_MR<- sQTL[!sQTL$SNP%in%df_GCST90002409$SNP,]
saveRDS(proxy_MR, 'missingFromOutcome.rds')
# Running the proxy search for SNP batch ---------------------------------------

token <- "2dd0e33b428c"


#proxy search
LDproxy_batch(proxy_MR$SNP, token = token, genome_build = "grch38", append = TRUE)