# Header ----------------------------------------------------------------------

# AUTHOR: John Eric Hamilton
# PURPOSE: Proxy Search for SNPs that are in the exposure and are missing the outcome
# NOTES: 
#   *Inputs:   - Outcome (df_GCST90002409)
#              - Proxy Data from ProxyBatch.R
#
#   *Outputs:  - outcomeToMerge with the processed proxied SNPs found (to be added to the outcome) 
            
library(TwoSampleMR)
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(reshape2)
library(LDlinkR)
# Loading Data and looking for missing SNPs in outcome ---------------------------------
### GWAS 'df_GCST90002409' (outcome) vs sQTL 'df' (exposure)
freadf <- function(df){return(as.data.frame(fread(df)))}

# Load datasets ----------------------------------------------------------------
df_GCST90002409 <- freadf("GCST90002409_buildGRCh38.tsv")

print(head(df_GCST90002409))

# Preping GWAS for LD Proxy ------------------------------------------------------
# Adding SNP column
df_GCST90002409$SNP <- paste("chr", df_GCST90002409$chromosome, ":", df_GCST90002409$base_pair_location, sep="")

df_GCST90002409$effect_allele <- toupper(df_GCST90002409$effect_allele)
df_GCST90002409$other_allele  <- toupper(df_GCST90002409$other_allele)
# Remove entry if there is multiple alleles in either the other_allele or effect_allele (indels)
df_GCST90002409 <- df_GCST90002409[df_GCST90002409$other_allele%in%c('A','T','C','G') & df_GCST90002409$effect_allele%in%c('A','T','C','G'),]
# Save
saveRDS(df_GCST90002409, 'VogelezangWithSNP.rds')
print(head(df_GCST90002409))
# Troubleshooting LD Proxy Errors ---------------------------------------------- 
proxy_MR <- readRDS('missingFromOutcome.rds')

print(head(proxy_MR))
# Running the proxy search -----------------------------------------------------
outcomeToMerge <- NULL
proxySearch <- read.delim('combined_query_snp_list_grch38.txt', row.names = NULL)

df_GCST90002409$SNP_Original <- df_GCST90002409$variant_id

token <- "2dd0e33b428c"
for(snp_p in proxy_MR$SNP){
  print(snp_p)
  #proxy search
  ldp=proxySearch[proxySearch$query_snp==snp_p,]
  #saveRDS(snp_p,"Data/Proxy/Proxy_GWAS_Vogelezang_snp.rds")
  if(sum(proxySearch$query_snp==snp_p)==0){
    next
  }
  
  ldp <- ldp[ldp$R2 > 0.8,]
  #ldp$Coord <- gsub("chr","",ldp$Coord)
  ldp=ldp[ldp$Coord%in%df_GCST90002409$SNP,]
  if(nrow(ldp)==0){
    next
    }
  ldp2=ldp[1,]
  toMerge=df_GCST90002409[df_GCST90002409$SNP==ldp2$Coord,]
  toMerge$SNP=ldp2$Coord
  toMerge$SNP_Original=ldp2$RS_Number
  if(ldp2$Distance!=0){
    ref=colsplit(ldp2$Correlated_Alleles,',',c('a','b'))
    alt=colsplit(ref$b,'=',c('exp','out'))
    alt[alt==TRUE]='T'
    ref=colsplit(ref$a,'=',c('exp','out'))
    ref[ref==TRUE]='T'
    # Check that both allele values exist before comparing
    if(length(toMerge$EA) == 0 || length(alt$out) == 0){
      cat("Allele information missing for SNP", snp_p, "\n")
      next
    }
    if(toMerge$EA==alt$out){
      toMerge$EA=alt$exp
      toMerge$Allele1=ref$exp
    }else{
      toMerge$EA=ref$exp
      toMerge$Allele1=alt$exp
    }
  }
  print("Binding...")
  outcomeToMerge=rbind(outcomeToMerge,toMerge)
  saveRDS(outcomeToMerge,"outcomeToMerge.rds")
}