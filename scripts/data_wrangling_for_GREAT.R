# Set working directory and load necessary packages
#setwd("C:/Users/rittere5/OneDrive - Michigan State University/Vitis-domatia/")
setwd("C:/Users/elean/OneDrive - Michigan State University/Vitis-domatia/")

######################## DATA WRANGLING FOR GREAT ########################

# Load data
dict <- read.csv("Vriparia-Vvinifera.tsv", sep='\t', header = FALSE)
df1 <- read.csv("DE_genes_Domatia_V_Leaf_588710.csv")
row.names(df1) <- df1$Gene
df2 <- read.csv("DE_genes_Domatia_V_Leaf_588711.csv")
row.names(df2) <- df2$Gene


# Shared between domatia and domatia
dom <- df1[rownames(df1)%in%rownames(df2),]

# Get Vvinifera orthologs to differentially expressed genes
vvi <- as.data.frame(merge(dict, dom, by.x = "V1", by.y = "Gene"))
#temp1 <- as.data.frame(merge(dict, dom, by.x = "V1", by.y = "Gene")) #Good for searching gene IDs
vvi <- vvi$V2
vvi <- as.data.frame(gsub(".t01", "", vvi))
colnames(vvi) <- c("V1")
vvi<-t(vvi)

# Copy data for GREAT
#writeClipboard(paste0(as.vector(vvi),collapse=","))

######################## DE ANALYSIS WITH GREAT DATA ########################
library(DESeq2)
library(dplyr)

# Load GREAT data into R and reformat
data1 <- read.csv(file = "GREAT_analysis/GREAT-data/GREAT_RPKM_samples2023-07-26_inflorescence.csv", row.names = "X")
data2 <- read.csv(file = "GREAT_analysis/GREAT-data/GREAT_RPKM_samples2023-07-26_leaves.csv", row.names = "X")
data3 <- read.csv(file = "GREAT_analysis/GREAT-data/GREAT_RPKM_samples2023-07-26_stem.csv", row.names = "X")
temp <- cbind(data1, data2, data3)
counts <- mutate_all(temp, function(x) as.numeric(as.integer(x)))
temp$Gene <- row.names(temp)

colData <- read.csv(file="GREAT_analysis/great_matrix.csv")
colData$tissue <- as.factor(colData$tissue)

##Making Col DESeq Data Set
dds = DESeqDataSetFromMatrix(countData = counts, colData = colData,  design = ~tissue)
#Make sure Control is the first level of treatment
dds$tissue <- relevel(dds$tissue, ref = "Control")
#Make sure that dds is properly formatted
as.data.frame(colData(dds))

##Run DESeq
dds <- DESeq(dds)
res <- results(dds)

## Double check result names
resultsNames(dds)

## Filter differentially expressed genes between inflorescences and leaves
deseq1 <- results(dds, contrast=c("tissue", "inflorescence", "Control"))
df.deseq1 <- as.data.frame(deseq1)
res1.deseq1<-df.deseq1[(df.deseq1$log2FoldChange>1|df.deseq1$log2FoldChange< -1)& !is.na(df.deseq1$log2FoldChange),]
inflor.leaves <- res1.deseq1[(res1.deseq1$padj<0.05)&!is.na(res1.deseq1$padj),]
inflor.leaves$Gene <- row.names(inflor.leaves)
inflor.upregulated <- inflor.leaves[(inflor.leaves$log2FoldChange>0),]

## Filter differentially expressed genes between inflorescences and stems
deseq1 <- results(dds, contrast=c("tissue", "inflorescence", "stem"))
df.deseq1 <- as.data.frame(deseq1)
res1.deseq1<-df.deseq1[(df.deseq1$log2FoldChange>1|df.deseq1$log2FoldChange< -1)& !is.na(df.deseq1$log2FoldChange),]
inflor.stem <- res1.deseq1[(res1.deseq1$padj<0.05)&!is.na(res1.deseq1$padj),]
inflor.stem$Gene <- row.names(inflor.stem)

## Filter differentially expressed genes between inflorescences and stems
deseq1 <- results(dds, contrast=c("tissue", "Control", "stem"))
df.deseq1 <- as.data.frame(deseq1)
res1.deseq1<-df.deseq1[(df.deseq1$log2FoldChange>1|df.deseq1$log2FoldChange< -1)& !is.na(df.deseq1$log2FoldChange),]
leaves.stem <- res1.deseq1[(res1.deseq1$padj<0.05)&!is.na(res1.deseq1$padj),]
leaves.stem$Gene <- row.names(leaves.stem)
stem.upregulated <- leaves.stem[(leaves.stem$log2FoldChange<0),]

# Remove any upregulated stem genes
test <- subset(inflor.upregulated, !(Gene %in% stem.upregulated$Gene))
