# Load the biomaRt library
library(biomaRt)

# Connect to the Ensembl BioMart database using the GRCh37 (hg19) archive.
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                dataset = "hsapiens_gene_ensembl", 
                host = "https://grch37.ensembl.org")

# Query the genomic attributes for the KNOP1 gene using its HGNC symbol.
knop1_location <- getBM(attributes = c("chromosome_name",
                                         "start_position",
                                         "end_position",
                                         "strand",
                                         "band"),
                        filters = "hgnc_symbol", 
                        values = "KNOP1", 
                        mart = mart)

# Display the results
print(knop1_location)
#build 38 KNOP1 16       19701937     19718235
#16:19,701,937-19,718,235
#build 37 Chromosome 16: 19,714,902-19,729,557
# plink --bfile Chr16.EUR --chr 16 --from-bp 19714902 --to-bp 19729557 --make-bed --out Chr16.EUR_subset
