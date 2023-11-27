# Set working directory and load necessary packages
setwd("C:/Users/rittere5/OneDrive - Michigan State University/Vitis-domatia/")
#setwd("C:/Users/elean/OneDrive - Michigan State University/Vitis-domatia/")

library("topGO")
library(dplyr)
library(tidyr)
library(ggplot2)

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

# For only upregulated genes
degenes2 <- degenes[degenes$log2FoldChange>0 ,]
final<-cbind(degenes2,all[match(degenes2$Gene,all$gene),])

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

top20.10a <- as.data.frame(top20.10[, c(1,3,6,8)])
top20.10a$Genotype <- c("588710")
colnames(top20.10a) <- c("GO.Term.Biological.Process", "DEGs", "Fold.Enrichment", "FDR.P.value", "Genotype")

top20.11a <- as.data.frame(top20.11[, c(1,3,6,8)])
top20.11a$Genotype <- c("588711")
colnames(top20.11a) <- c("GO.Term.Biological.Process", "DEGs", "Fold.Enrichment", "FDR.P.value", "Genotype")

df10.a <- as.data.frame(df10[df10$upload_1..over.under..588710=="+" ,] [, c(1,3,6,8)])
df10.a$Genotype <- c("588710")
colnames(df10.a) <- c("GO.Term.Biological.Process", "DEGs", "Fold.Enrichment", "FDR.P.value", "Genotype")

df11.a <- as.data.frame(df11[df11$upload_1..over.under..588711=="+" ,] [, c(1,3,6,8)])
df11.a$Genotype <- c("588711")
colnames(df11.a) <- c("GO.Term.Biological.Process", "DEGs", "Fold.Enrichment", "FDR.P.value", "Genotype")

overlap1 <- df10.a[df10.a$GO.Term.Biological.Process %in% top20.11a$GO.Term.Biological.Process ,]
overlap2 <- df11.a[df11.a$GO.Term.Biological.Process %in% top20.10a$GO.Term.Biological.Process ,]

data <- unique(rbind(overlap1, overlap2, top20.10a, top20.11a))

ggplot(data, (aes(x=Genotype, y=GO.Term.Biological.Process, color = as.numeric(FDR.P.value), size=as.numeric(DEGs)))) + 
  geom_point(shape = 19) +
  scale_color_gradient(low = "orange", high = "blue") +
  theme_classic() +
  scale_y_discrete(limits=rev) +
  ylab("GO Term") 

# Look at overlap for all significant and positively enriched GO terms
sign10 <- df10.a[df10.a$FDR.P.value <0.05 ,]
sign11 <- df11.a[df11.a$FDR.P.value<0.05 ,]
tempa <- sign10[sign10$GO.Term.Biological.Process %in% sign11$GO.Term.Biological.Process ,]
tempb <- sign11[sign11$GO.Term.Biological.Process %in% sign10$GO.Term.Biological.Process ,]
sign.overlap <- unique(rbind(tempa, tempb))

ggplot(sign.overlap, (aes(x=Genotype, y=GO.Term.Biological.Process, color = as.numeric(FDR.P.value), size=as.numeric(Fold.Enrichment)))) + 
  geom_point() +
  scale_color_gradient(low = "orange", high = "blue") +
  theme_classic()

# # With top fold changes
# top20.10 <- head(sign10[order(sign10$Fold.Enrichment, decreasing = TRUE),], n = 20)
# top20.11 <- head(sign11[order(sign11$Fold.Enrichment, decreasing = TRUE),], n = 20)
# 
# df10.a <- as.data.frame(df10[df10$upload_1..over.under..588710=="+" ,] [, c(1,6,8)])
# df10.a$Genotype <- c("588710")
# colnames(df10.a) <- c("GO.Term.Biological.Process", "Fold.Enrichment", "FDR.P.value", "Genotype")
# 
# df11.a <- as.data.frame(df11[df11$upload_1..over.under..588711=="+" ,] [, c(1,6,8)])
# df11.a$Genotype <- c("588711")
# colnames(df11.a) <- c("GO.Term.Biological.Process", "Fold.Enrichment", "FDR.P.value", "Genotype")
# 
# overlap1 <- df10.a[df10.a$GO.Term.Biological.Process %in% top20.11$GO.Term.Biological.Process ,]
# overlap2 <- df11.a[df11.a$GO.Term.Biological.Process %in% top20.10$GO.Term.Biological.Process ,]
# 
# data <- unique(rbind(overlap1, overlap2, top20.10, top20.11))
# 
# ggplot(data, (aes(x=Genotype, y=GO.Term.Biological.Process, color = as.numeric(FDR.P.value), size=as.numeric(Fold.Enrichment)))) + 
#   geom_point() +
#   scale_color_gradient(low = "orange", high = "blue") +
#   theme_classic()

######################## MAKE PLOT FROM PANTHER RESULTS (W BONFERRONI) FOR BOTH GENOTYPES ########################
df10 <- read.csv("GO-term-enrichment/DE_genes_Domatia_V_Leaf_588710_GO_term_enrichment_analysis_bonferroni.csv", skip = 12, sep = '\t')
colnames(df10) <- paste( colnames(df10), ".588710", sep="")
df11 <- read.csv("GO-term-enrichment/DE_genes_Domatia_V_Leaf_588711_GO_term_enrichment_analysis_bonferroni.csv", skip = 12, sep = '\t')
colnames(df11) <- paste( colnames(df11), ".588711", sep="")

# top20.10 <- head(df10[order(df10[df10$upload_1..over.under..588710=="+",]$upload_1..FDR..588710),], n = 20)
# top20.11 <- head(df11[order(df11[df11$upload_1..over.under..588711=="+",]$upload_1..FDR..588711),], n = 20)
# 
# top20.10a <- as.data.frame(top20.10[, c(1,3,6,8)])
# top20.10a$Genotype <- c("588710")
# colnames(top20.10a) <- c("GO.Term.Biological.Process", "DEGs", "Fold.Enrichment", "FDR.P.value", "Genotype")
# 
# top20.11a <- as.data.frame(top20.11[, c(1,3,6,8)])
# top20.11a$Genotype <- c("588711")
# colnames(top20.11a) <- c("GO.Term.Biological.Process", "DEGs", "Fold.Enrichment", "FDR.P.value", "Genotype")

df10.a <- as.data.frame(df10[df10$upload_1..over.under..588710=="+" ,] [, c(1,3,5,6,7)])
df10.a$Genotype <- c("588710")
colnames(df10.a) <- c("GO.Term.Biological.Process", "DEGs", "Over.Under", "Fold.Enrichment", "P.value", "Genotype")

df11.a <- as.data.frame(df11[df11$upload_1..over.under..588711=="+" ,] [, c(1,3,5,6,7)])
df11.a$Genotype <- c("588711")
colnames(df11.a) <- c("GO.Term.Biological.Process", "DEGs", "Over.Under", "Fold.Enrichment", "P.value", "Genotype")

# overlap1 <- df10.a[df10.a$GO.Term.Biological.Process %in% top20.11a$GO.Term.Biological.Process ,]
# overlap2 <- df11.a[df11.a$GO.Term.Biological.Process %in% top20.10a$GO.Term.Biological.Process ,]
# 
# data <- unique(rbind(overlap1, overlap2, top20.10a, top20.11a))
# 
# ggplot(data, (aes(x=Genotype, y=GO.Term.Biological.Process, color = as.numeric(FDR.P.value), size=as.numeric(DEGs)))) + 
#   geom_point(shape = 19) +
#   scale_color_gradient(low = "orange", high = "blue") +
#   theme_classic() +
#   scale_y_discrete(limits=rev) +
#   ylab("GO Term") 

# Look at overlap for all significant and positively enriched GO terms
sign10 <- df10.a[df10.a$P.value <0.05 ,]
sign11 <- df11.a[df11.a$P.value<0.05 ,]
tempa <- sign10[sign10$GO.Term.Biological.Process %in% sign11$GO.Term.Biological.Process ,]
tempb <- sign11[sign11$GO.Term.Biological.Process %in% sign10$GO.Term.Biological.Process ,]
sign.overlap <- unique(rbind(tempa, tempb))

ggplot(sign.overlap, (aes(x=Genotype, y=GO.Term.Biological.Process, color = as.numeric(Fold.Enrichment), size=as.numeric(DEGs)))) + 
  geom_point() +
  scale_color_gradient(low = "orange", high = "blue", name = "Fold Enrichment") +
  scale_y_discrete(limits=rev) +
  theme_classic() +
  ylab("GO Term (Biological Process)") +
  scale_size_continuous(name = "DEGs") +
  theme(axis.text.x = element_text(size=11, face="bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=11, face="bold"))

######################## MAKE PLOT FROM PANTHER RESULTS (W BONFERRONI) FOR BOTH GENOTYPES - CC ########################
df10 <- read.csv("GO-term-enrichment/DE_genes_Domatia_V_Leaf_588710_GO_term_enrichment_analysis_bonferroni_CC.csv", skip = 12, sep = '\t')
colnames(df10) <- paste( colnames(df10), ".588710", sep="")
df11 <- read.csv("GO-term-enrichment/DE_genes_Domatia_V_Leaf_588711_GO_term_enrichment_analysis_bonferroni_CC.csv", skip = 12, sep = '\t')
colnames(df11) <- paste( colnames(df11), ".588711", sep="")

df10.a <- as.data.frame(df10[df10$upload_1..over.under..588710=="+" ,] [, c(1,3,5,6,7)])
df10.a$Genotype <- c("588710")
colnames(df10.a) <- c("GO.Term.Cellular.Component", "DEGs", "Over.Under", "Fold.Enrichment", "P.value", "Genotype")

df11.a <- as.data.frame(df11[df11$upload_1..over.under..588711=="+" ,] [, c(1,3,5,6,7)])
df11.a$Genotype <- c("588711")
colnames(df11.a) <- c("GO.Term.Cellular.Component", "DEGs", "Over.Under", "Fold.Enrichment", "P.value", "Genotype")

# Look at overlap for all significant and positively enriched GO terms
sign10 <- df10.a[df10.a$P.value <0.05 ,]
sign11 <- df11.a[df11.a$P.value<0.05 ,]
tempa <- sign10[sign10$GO.Term.Cellular.Component %in% sign11$GO.Term.Cellular.Component ,]
tempb <- sign11[sign11$GO.Term.Cellular.Component %in% sign10$GO.Term.Cellular.Component ,]
sign.overlap <- unique(rbind(tempa, tempb))

sign.overlap1 <- sign.overlap[!grepl("Unclassified*", sign.overlap$GO.Term.Cellular.Component),]

ggplot(sign.overlap1, (aes(x=Genotype, y=GO.Term.Cellular.Component, color = as.numeric(Fold.Enrichment), size=as.numeric(DEGs)))) +
  geom_point() +
  scale_color_gradient(low = "orange", high = "blue", name = "Fold Enrichment") +
  scale_y_discrete(limits=rev) +
  theme_classic() +
  ylab("GO Term (Cellular Component)") +
  scale_size_continuous(name = "DEGs") +
  theme(axis.text.x = element_text(size=11, face="bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=11, face="bold"))

# #Plotting all significant terms (even if they DON'T overlap!)
# overlap1 <- df10.a[df10.a$GO.Term.Cellular.Component %in% sign11$GO.Term.Cellular.Component ,]
# overlap2 <- df11.a[df11.a$GO.Term.Cellular.Component %in% sign10$GO.Term.Cellular.Component ,]
# 
# data <- unique(rbind(overlap1, overlap2, sign10, sign11))
# 
# ggplot(data, (aes(x=Genotype, y=GO.Term.Cellular.Component, color = as.numeric(P.value), size=as.numeric(DEGs)))) +
#   geom_point(shape = 19) +
#   scale_color_gradient(low = "orange", high = "blue") +
#   theme_classic() +
#   scale_y_discrete(limits=rev) +
#   ylab("GO Term (Cellular Component)")


######################## INVESTIGATING GO TERM ENRICHMENT FOR CONTROL V CONTROL - BP ########################

data.bp <- read.csv("GO-term-enrichment/DE_genes_Control_588710_V_Control_588711_GO_term_enrichment_analysis_bonferroni.csv", skip = 12, sep = '\t')
sign.bp <- data.bp[data.bp$upload_1..P.value. <0.05 ,]
sign.bp.over <- sign.bp[sign.bp$upload_1..over.under.=="+" ,]

# Notes: a good amount of overlap with what is different between domatia and control - interesting! Quantified below

# df10 <- read.csv("GO-term-enrichment/DE_genes_Domatia_V_Leaf_588710_GO_term_enrichment_analysis.csv", skip = 11, sep = '\t')
# colnames(df10) <- paste( colnames(df10), ".588710", sep="")
# df11 <- read.csv("GO-term-enrichment/DE_genes_Domatia_V_Leaf_588711_GO_term_enrichment_analysis.csv", skip = 11, sep = '\t')
# colnames(df11) <- paste( colnames(df11), ".588711", sep="")
# 
# top20.10 <- head(df10[order(df10[df10$upload_1..over.under..588710=="+",]$upload_1..FDR..588710),], n = 20)
# top20.11 <- head(df11[order(df11[df11$upload_1..over.under..588711=="+",]$upload_1..FDR..588711),], n = 20)
# 
# top20.10a <- as.data.frame(top20.10[, c(1,3,6,8)])
# top20.10a$Genotype <- c("588710")
# colnames(top20.10a) <- c("GO.Term.Biological.Process", "DEGs", "Fold.Enrichment", "FDR.P.value", "Genotype")
# 
# top20.11a <- as.data.frame(top20.11[, c(1,3,6,8)])
# top20.11a$Genotype <- c("588711")
# colnames(top20.11a) <- c("GO.Term.Biological.Process", "DEGs", "Fold.Enrichment", "FDR.P.value", "Genotype")
# 
# df10.a <- as.data.frame(df10[df10$upload_1..over.under..588710=="+" ,] [, c(1,3,6,8)])
# df10.a$Genotype <- c("588710")
# colnames(df10.a) <- c("GO.Term.Biological.Process", "DEGs", "Fold.Enrichment", "FDR.P.value", "Genotype")
# 
# df11.a <- as.data.frame(df11[df11$upload_1..over.under..588711=="+" ,] [, c(1,3,6,8)])
# df11.a$Genotype <- c("588711")
# colnames(df11.a) <- c("GO.Term.Biological.Process", "DEGs", "Fold.Enrichment", "FDR.P.value", "Genotype")
# 
# # Look at overlap for all significant and positively enriched GO terms
# sign10 <- df10.a[df10.a$FDR.P.value <0.05 ,]
# sign11 <- df11.a[df11.a$FDR.P.value<0.05 ,]
# tempa <- sign10[sign10$GO.Term.Biological.Process %in% sign11$GO.Term.Biological.Process ,]
# tempb <- sign11[sign11$GO.Term.Biological.Process %in% sign10$GO.Term.Biological.Process ,]
# sign.overlap <- unique(rbind(tempa, tempb))
# 
# bp.overlap <- sign.bp.over[sign.bp.over$GO.biological.process.complete %in% sign.overlap$GO.Term.Biological.Process ,]

# # 32 of the GO terms overlap with those shared between domatia v control datasets

######################## INVESTIGATING GO TERM ENRICHMENT FOR CONTROL V CONTROL - CC ########################

data.cc <- read.csv("GO-term-enrichment/DE_genes_Control_588710_V_Control_588711_GO_term_enrichment_analysis_bonferroni_CC.csv", skip = 12, sep = '\t')
sign.cc <- data.cc[data.cc$upload_1..P.value. <0.05 ,]
sign.cc.over <- sign.cc[sign.cc$upload_1..over.under.=="+" ,]

# # Same as above - let's quantify overlap with those shared and significant in domatia v control datasets
# df10 <- read.csv("GO-term-enrichment/DE_genes_Domatia_V_Leaf_588710_GO_term_enrichment_analysis_bonferroni_CC.csv", skip = 12, sep = '\t')
# colnames(df10) <- paste( colnames(df10), ".588710", sep="")
# df11 <- read.csv("GO-term-enrichment/DE_genes_Domatia_V_Leaf_588711_GO_term_enrichment_analysis_bonferroni_CC.csv", skip = 12, sep = '\t')
# colnames(df11) <- paste( colnames(df11), ".588711", sep="")
# 
# df10.a <- as.data.frame(df10[df10$upload_1..over.under..588710=="+" ,] [, c(1,3,5,6,7)])
# df10.a$Genotype <- c("588710")
# colnames(df10.a) <- c("GO.Term.Cellular.Component", "DEGs", "Over.Under", "Fold.Enrichment", "P.value", "Genotype")
# 
# df11.a <- as.data.frame(df11[df11$upload_1..over.under..588711=="+" ,] [, c(1,3,5,6,7)])
# df11.a$Genotype <- c("588711")
# colnames(df11.a) <- c("GO.Term.Cellular.Component", "DEGs", "Over.Under", "Fold.Enrichment", "P.value", "Genotype")
# 
# # Look at overlap for all significant and positively enriched GO terms
# sign10 <- df10.a[df10.a$P.value <0.05 ,]
# sign11 <- df11.a[df11.a$P.value<0.05 ,]
# tempa <- sign10[sign10$GO.Term.Cellular.Component %in% sign11$GO.Term.Cellular.Component ,]
# tempb <- sign11[sign11$GO.Term.Cellular.Component %in% sign10$GO.Term.Cellular.Component ,]
# sign.overlap <- unique(rbind(tempa, tempb))
# 
# cc.overlap <- sign.cc.over[sign.cc.over$GO.cellular.component.complete %in% sign.overlap$GO.Term.Cellular.Component ,]
# 
# # 6 of the GO terms overlap with those shared between domatia v control datasets - all cell wall related

######################## INVESTIGATING CELL WALL GENES FROM CONTROL V CONTROL ########################
cell.wall <- read.csv("GO-term-enrichment/DE_genes_Control_588710_V_Control_588711_plant-type_cell_wall_gene_list.csv", sep = "\t", header = FALSE)
data10 <- read.csv("differentially_expressed_genes/DE_genes_Domatia_V_Leaf_588710.csv")
data11 <- read.csv("differentially_expressed_genes/DE_genes_Domatia_V_Leaf_588711.csv")

overlap10 <- merge(cell.wall, data10, by.x="V2", by.y="ensembl_gene_id")
overlap11 <- merge(cell.wall, data11, by.x="V2", by.y="ensembl_gene_id")

both <- merge(overlap10, overlap11, by="V2")
