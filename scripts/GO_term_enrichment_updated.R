# Set working directory and load necessary packages
#setwd("C:/Users/rittere5/OneDrive - Michigan State University/Vitis-domatia/")
setwd("C:/Users/elean/OneDrive - Michigan State University/Vitis-domatia/")

library("topGO")
library(dplyr)
library(tidyr)

######################## GO TERM ENRICHMENT WITH DOMATIA V LEAF IN 588710 ######################## 

# Load data

## Load functional annotations
all <- read.csv("Vriparia-functional-annotations.tsv", sep='\t')

## Add gene data to all if needed
dict <- read.csv("cds-to-gene.tsv", sep = '\t')
colnames(dict) <- c("Transcript", "gene")
all <- merge(all, dict, by="Transcript", all = TRUE)
all <- all[,c(1, 8, 2:7)]
all <- all %>% drop_na(gene)
all <- all[!duplicated(all),]

## Load differentially expressed genes
degenes <- read.csv("differentially_expressed_genes/DE_genes_Domatia_V_Leaf_588710.csv")

# Modify DE gene list
final<-cbind(degenes,all[match(degenes$Gene,all$gene),])

# Copy Arabidopsis orthologs
writeClipboard(paste0(as.vector(gsub("\\.1", "", final$Arabidopsis_blast_hit))))

               