# Set working directory and load necessary packages
setwd("C:/Users/rittere5/OneDrive - Michigan State University/Vitis-domatia/")
#setwd("C:/Users/elean/OneDrive - Michigan State University/Vitis-domatia/")

library(dplyr)
library(tidyr)
library(ggplot2)

# Read in differentially expressed genes
dom710 <- read.csv("differentially_expressed_genes/DE_genes_Domatia_V_Leaf_588710.csv")
dom711 <- read.csv("differentially_expressed_genes/DE_genes_Domatia_V_Leaf_588711.csv")

# Read in EFN transport genes
efn <- read.csv("efn-genes.csv")

# Get overlap for each genotype individually
shared710 <- merge(efn, dom710, by.x="A.thaliana.ID",by.y="ensembl_gene_id")
shared711 <- merge(efn, dom711, by.x="A.thaliana.ID",by.y="ensembl_gene_id")
