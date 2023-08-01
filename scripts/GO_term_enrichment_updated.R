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
writeClipboard(paste0(as.vector(na.omit(gsub("\\.[0-9]", "", final$Arabidopsis_blast_hit)))))

rm(list = ls())

######################## GO TERM ENRICHMENT WITH DOMATIA V LEAF IN 588711 ######################## 

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
degenes <- read.csv("differentially_expressed_genes/DE_genes_Domatia_V_Leaf_588711.csv")

# Modify DE gene list
final<-cbind(degenes,all[match(degenes$Gene,all$gene),])

# Copy Arabidopsis orthologs
writeClipboard(paste0(as.vector(na.omit(gsub("\\.[0-9]", "", final$Arabidopsis_blast_hit)))))  

rm(list = ls())

######################## GO TERM ENRICHMENT WITH DOMATIA V DOMATIA ######################## 

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
degenes <- read.csv("differentially_expressed_genes/DE_genes_Domatia_710_v_711.csv")

# Modify DE gene list
final<-cbind(degenes,all[match(degenes$Gene,all$gene),])

# Copy Arabidopsis orthologs
writeClipboard(paste0(as.vector(na.omit(gsub("\\.[0-9]", "", final$Arabidopsis_blast_hit)))))               

rm(list = ls())

######################## GO TERM ENRICHMENT WITH CONTROL V CONTROL ######################## 

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
degenes <- read.csv("differentially_expressed_genes/DE_genes_Control_710_v_711.csv")

# Modify DE gene list
final<-cbind(degenes,all[match(degenes$Gene,all$gene),])

# Copy Arabidopsis orthologs
writeClipboard(paste0(as.vector(na.omit(gsub("\\.[0-9]", "", final$Arabidopsis_blast_hit)))))               

rm(list = ls())
######################## MAKE PLOT FROM PANTHER RESULTS FOR BOTH GENOTYPES ########################
df10 <- read.csv("GO-term-enrichment/DE_genes_Domatia_V_Leaf_588710_GO_term_enrichment_analysis.csv", skip = 11, sep = '\t')
colnames(df10) <- paste( colnames(df10), ".588710", sep="")
df11 <- read.csv("GO-term-enrichment/DE_genes_Domatia_V_Leaf_588711_GO_term_enrichment_analysis.csv", skip = 11, sep = '\t')
colnames(df11) <- paste( colnames(df11), ".588711", sep="")

top20.10 <- head(df10[order(df10[df10$upload_1..over.under..588710=="+",]$upload_1..FDR..588710),], n = 20)
top20.11 <- head(df11[order(df11[df11$upload_1..over.under..588711=="+",]$upload_1..FDR..588711),], n = 20)

overlap1 <- merge(top20.10, df11, by.x="GO.biological.process.complete.588710", by.y="GO.biological.process.complete.588711")
overlap2 <- merge(df10, top20.11, by.x="GO.biological.process.complete.588710", by.y="GO.biological.process.complete.588711")

data <- unique(rbind(overlap1, overlap2))
data10 <- data[, 1:8]
data11 <- data[, c(1,9:15)]

final <- data[, c(1,6,8,13,15)]
colnames(final) <- c("GO Term - Biological Process", "588710 Fold Change", "588710 FDR P-value",
                     "588711 Fold Change", "588711 FDR P-value")

sign10 <- df10[df10$upload_1..FDR..588710<0.05 ,]
sign11 <- df11[df11$upload_1..FDR..588711<0.05 ,]
sign.overlap <-merge(sign10, sign11, by.x="GO.biological.process.complete.588710", by.y="GO.biological.process.complete.588711")
