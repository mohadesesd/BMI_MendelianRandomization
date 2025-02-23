# Header ----------------------------------------------------------------------

# Adapted from John Eric Hamilton, modified by Mohadese Sayahian Dehkordi
# PURPOSE: Plotting results from the Data collection 
# NOTES: 
#   *Inputs:   - MR Results filtered and non-filtered 
#
#   *Outputs:  - Forest Plot of significant genes based on OR 
#              - Gene Name Conversion Table (from Gene Ensemble ID to common gene name)

# Library and Data Loading ----------------------------------------------------
library(tidyverse)
library(stringr)
# Load the filtered (for significant p-value) MR results
result_signif_total <- readRDS('./Pancreas_MR/allResultSignif.rds')
filtered <- readRDS('./Pancreas_MR/allResults_filtered.rds')
filtered_MHC <- readRDS('./Pancreas_MR/allResults_filtered_MHC_region.rds')

# Forest Plot statistics -------------------------------------------------------

# Calculating the Odd Ratio
result_signif_total$OR <- exp(result_signif_total$b)

# Calculating the confidence interval
result_signif_total$CI_lower <- result_signif_total$b -  result_signif_total$se*qnorm(0.975)
result_signif_total$CI_upper <- result_signif_total$b +  result_signif_total$se*qnorm(0.975)

# Order properly and Create factors for ordering the entries for ggplot
result_signif_total <- result_signif_total[order(result_signif_total$b),]

saveRDS(result_signif_total, './Pancreas_MR/allResults_OR.rds')


### Filtering result_signif_total to get only non-MHC region
filtered <- result_signif_total[(result_signif_total$id.exposure %in% filtered$id.exposure),]

filtered <- filtered[order(filtered$b),]

saveRDS(filtered, './Pancreas_MR/allResults_filteredWithNames.rds')

### Filtering result_signif_total to get only MHC region
filtered_MHC <- result_signif_total[(result_signif_total$id.exposure %in% filtered_MHC$id.exposure),]

filtered_MHC <- filtered_MHC[order(filtered_MHC$b),]

saveRDS(filtered_MHC, './Pancreas_MR/allResults_filteredMHCWithNames.rds')


# Remove bad columns (duplicated geneName.1 for some reason) - (Might not always need)
#result_signif_total <- subset(result_signif_total, select=-c(15))
#filtered_WithNames <- subset(filtered_WithNames, select=-c(15))
#filtered_MHC_WithNames <- subset(filtered_MHC_WithNames, select=-c(15))

gtf_reduced <- readRDS('GeneCode_ConversionTable.rds')


# Creating a column in the result_signif_total table to make the conversion ----
result_signif_total_WithNames <- data.frame()

for (j in 1:dim(result_signif_total)[1]){
  #entry <- gtf_reduced %>% filter(gtf_reduced$id == strsplit(result_signif_total$id.exposure[j], split='\\.')[[1]][1])
  entry <- gtf_reduced[gtf_reduced$id == sub(".*:(ENSG[0-9]+)\\..*", "\\1", result_signif_total$id.exposure[j]),]
  name <- entry$name
  
  # To fix bug when the ensembleId is not in the GTF Conversion table
  if (dim(entry)[1] == 0){
    print("Not found")
    name <- result_signif_total$id.exposure[j]
  }
  
  temp_name <- cbind(result_signif_total[j,], name)
  result_signif_total_WithNames <- rbind(result_signif_total_WithNames, temp_name) 
  print("found")
}

# Change column Name
names(result_signif_total_WithNames)[names(result_signif_total_WithNames) == "name"] <- "geneName"
result_signif_total <- result_signif_total_WithNames
print(head(result_signif_total))
saveRDS(result_signif_total, './Pancreas_MR/allResultSignifWithNames.rds')


# Making a forest plot ---------------------------------------------------------
# 3 plots for total_signif, MHC, and non-MHC plots

data <- result_signif_total
###ALL
p_all <- ggplot(data, aes(x = b, y = geneName, xmin = CI_lower, xmax = CI_upper)) +
  geom_vline(xintercept = 0) +
  geom_point() +
  geom_errorbarh() + 
  xlab("Effect") + 
  ylab("Gene Name") +
  ggtitle("All Genes")

# Save Plot
pdf('./Pancreas_MR/AllSignifForestPlot.pdf', height = 12, width = 10)

#Show Plot
p_all
dev.off()

# Creating a column in the result_signif_total table to make the conversion ----
filtered_MHC_WithNames <- data.frame()

for (j in 1:dim(filtered_MHC)[1]){
  #entry <- gtf_reduced %>% filter(gtf_reduced$id == strsplit(result_signif_total$id.exposure[j], split='\\.')[[1]][1])
  entry <- gtf_reduced[gtf_reduced$id == sub(".*:(ENSG[0-9]+)\\..*", "\\1", filtered_MHC$id.exposure[j]),]
  name <- entry$name
  
  # To fix bug when the ensembleId is not in the GTF Conversion table
  if (dim(entry)[1] == 0){
    print("Not found")
    name <- filtered_MHC$id.exposure[j]
  }
  
  temp_name <- cbind(filtered_MHC[j,], name)
  filtered_MHC_WithNames <- rbind(filtered_MHC_WithNames, temp_name) 
  print("found")
}

# Change column Name
names(filtered_MHC_WithNames)[names(filtered_MHC_WithNames) == "name"] <- "geneName"
filtered_MHC <- filtered_MHC_WithNames
print(head(filtered_MHC))


###MHC
p_MHC <- ggplot(filtered_MHC_WithNames, aes(x = b, y = geneName, xmin = CI_lower, xmax = CI_upper)) +
  geom_vline(xintercept = 0) +
  geom_point() +
  geom_errorbarh() + 
  xlab("Effect") + 
  ylab("Gene Name") +
  ggtitle("All Genes MHC Region")

# Save Plot
pdf('./Pancreas_MR/AllSignif_MHCForestPlot.pdf', height = 12, width = 10)

#Show Plot
p_MHC
dev.off()


# Creating a column in the result_signif_total table to make the conversion ----
filtered_WithNames <- data.frame()

for (j in 1:dim(filtered)[1]){
  #entry <- gtf_reduced %>% filter(gtf_reduced$id == strsplit(result_signif_total$id.exposure[j], split='\\.')[[1]][1])
  entry <- gtf_reduced[gtf_reduced$id == sub(".*:(ENSG[0-9]+)\\..*", "\\1", filtered$id.exposure[j]),]
  name <- entry$name
  
  # To fix bug when the ensembleId is not in the GTF Conversion table
  if (dim(entry)[1] == 0){
    print("Not found")
    name <- filtered$id.exposure[j]
  }
  
  temp_name <- cbind(filtered[j,], name)
  filtered_WithNames <- rbind(filtered_WithNames, temp_name) 
  print("found")
}

# Change column Name
names(filtered_WithNames)[names(filtered_WithNames) == "name"] <- "geneName"
filtered <- filtered_WithNames
print(head(filtered))

### NON-MHC
p_nonMHC <- ggplot(filtered_WithNames, aes(x = b, y = geneName, xmin = CI_lower, xmax = CI_upper)) +
  geom_vline(xintercept = 0) +
  geom_point() +
  geom_errorbarh() + 
  xlab("Effect") + 
  ylab("Gene Name") +
  ggtitle("All Genes Non MHC") 

# Save Plot
pdf('./Pancreas_MR/AllSignif_nonMHCForestPlot.pdf', height = 12, width = 10)

#Show Plot
p_nonMHC
dev.off()
