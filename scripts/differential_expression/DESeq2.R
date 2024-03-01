# Set working directory and load necessary packages
#setwd("C:/Users/rittere5/OneDrive - Michigan State University/Vitis-domatia/")
setwd("C:/Users/elean/OneDrive - Michigan State University/Vitis-domatia/")
library(DESeq2)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(ggh4x)
library(reshape2)
library(gplots)
library(RColorBrewer)

############################################# RUNNING DESEQ2 #############################################
###Preparing Col files for DESeq 
ff <- list.files( path = "./star", pattern = "*ReadsPerGene.out.tab$", full.names = TRUE )
counts.files <- lapply( ff, read.table, skip = 4 )
counts <- as.data.frame( sapply( counts.files, function(x) x[ , 4 ] ) )
ff <- gsub( "ReadsPerGene[.]out[.]tab", "", ff )
ff <- gsub( "[.]/star/", "", ff )
colnames(counts) <- ff
row.names(counts) <- counts.files[[1]]$V1
colData <- read.csv(file="matrix.csv")
colData$genotype <- as.factor(colData$genotype)

##Making Col DESeq Data Set
dds = DESeqDataSetFromMatrix(countData = counts, colData = colData,  design = ~genotype + treatment + genotype:treatment)
#Make sure Control is the first level of treatment
dds$treatment <- relevel(dds$treatment, ref = "Control")
#Make sure that dds is properly formatted
as.data.frame(colData(dds))

##Run DESeq
dds <- DESeq(dds)
res <- results(dds)

############################### EXTRACTING AND FILTERING SPECIFIC SAMPLES ################################
## Double check result names
resultsNames(dds)

## Filter differentially expressed genes between 588710 Domatia and Leaf
deseq1 <- results(dds, contrast=c("treatment", "Domatia", "Control"))
df.deseq1 <- as.data.frame(deseq1)
res1.deseq1<-df.deseq1[(df.deseq1$log2FoldChange>1|df.deseq1$log2FoldChange< -1)& !is.na(df.deseq1$log2FoldChange),]
DOM.710 <- res1.deseq1[(res1.deseq1$padj<0.05)&!is.na(res1.deseq1$padj),]

## Filter differentially expressed genes between 588711 Domatia and Leaf
deseq2 <- results(dds, list(c("treatment_Domatia_vs_Control","genotype588711.treatmentDomatia")))
df.deseq2 <- as.data.frame(deseq2)
res1.deseq2<-df.deseq2[(df.deseq2$log2FoldChange>1|df.deseq2$log2FoldChange< -1)& !is.na(df.deseq2$log2FoldChange),]
DOM.711 <- res1.deseq2[(res1.deseq2$padj<0.05)&!is.na(res1.deseq2$padj),]

## What genes are different between domatia and leaf between the two genotypes?
deseq3 <- results(dds, name="genotype588711.treatmentDomatia")
df.deseq3 <- as.data.frame(deseq3)
res1.deseq3<-df.deseq3[(df.deseq3$log2FoldChange>1|df.deseq3$log2FoldChange< -1)& !is.na(df.deseq3$log2FoldChange),]
DOM.geno <- res1.deseq3[(res1.deseq3$padj<0.05)&!is.na(res1.deseq3$padj),]

## What genes are different between genotypes for the control leaf tissue?
genores <- results(dds, name="genotype_588711_vs_588710")
geno <- as.data.frame(genores)
geno1<-geno[(geno$log2FoldChange>1|geno$log2FoldChange< -1)& !is.na(geno$log2FoldChange),]
geno2 <- geno1[(geno1$padj<0.05)&!is.na(geno1$padj),]

############################### ATHA ORTHOLOGS ###############################
# Load diamond tsv file for Vitis riparia and Arabidopsis thaliana
vriatha <- read.csv("Vriparia-Athaliana.tsv", header = FALSE, sep = "\t")
vriatha1 <- vriatha %>% distinct(V1, .keep_all = TRUE)
row.names(vriatha1) <- vriatha1$V1

### DOM.GENO
# Get A. thaliana orthologs for differentially expressed genes
match <- merge(vriatha1,DOM.geno, by=0, all="T")
match2 <- match[!is.na(match$baseMean),]
row.names(match2) <- match2$Row.names
match3 <- match2[c(2, 3, 14, 15, 16, 17, 18, 19)]
match3$Gene <- row.names(match3)
matchF <- match3[c(9, 2, 3, 4, 5, 6, 7, 8)]
names(matchF)[names(matchF) == 'V2'] <- 'Athaliana.ortholog'

# Modify file for biomaRt
col2 <- matchF$Athaliana.ortholog
test <- sapply(col2,function(ii) gsub("\\..*","",ii))
matchF$Athaliana.ortholog <- test
rm(col2,test)

library(biomaRt)

# Read in genes
atha <- matchF$Athaliana.ortholog
atha <- paste(atha, ".1", sep="")

# Load Ensembl Plants into biomaRt
ensembl_plants = useMart(biomart="plants_mart", host = "plants.ensembl.org")

# Select dataset to use
datasets <- listDatasets(ensembl_plants)
ensembl = useDataset("athaliana_eg_gene", ensembl_plants)

# Figure out what filters and attributes will be used
#filters = listFilters(ensembl)
#filters[25,]
#attributes = listAttributes(ensembl)
#attributes[5,]

# Use getBM to get attributes for values
results1 <- getBM(attributes=c('ensembl_gene_id', 'description','tair_symbol', 'name_1006', 'definition_1006'),
                  filters = 'tair_locus_model',           
                  values = atha, 
                  mart = ensembl)

# Combine with vri genes
DOM.geno.ann <- cbind(matchF,results1[match(matchF$Athaliana.ortholog,results1$ensembl_gene_id),])
rm(datasets, ensembl, ensembl_plants, results1, atha)

### DOM.710
# Get A. thaliana orthologs for differentially expressed genes in domatia between genotypes
match <- merge(vriatha1,DOM.710, by=0, all="T")
match2 <- match[!is.na(match$baseMean),]
row.names(match2) <- match2$Row.names
match3 <- match2[c(2, 3, 14, 15, 16, 17, 18, 19)]
match3$Gene <- row.names(match3)
matchF <- match3[c(9, 2, 3, 4, 5, 6, 7, 8)]
names(matchF)[names(matchF) == 'V2'] <- 'Athaliana.ortholog'

# Modify file for biomaRt
col2 <- matchF$Athaliana.ortholog
test <- sapply(col2,function(ii) gsub("\\..*","",ii))
matchF$Athaliana.ortholog <- test
rm(col2,test)

library(biomaRt)

# Read in genes
atha <- matchF$Athaliana.ortholog
atha <- paste(atha, ".1", sep="")

# Load Ensembl Plants into biomaRt
ensembl_plants = useMart(biomart="plants_mart", host = "plants.ensembl.org")

# Select dataset to use
datasets <- listDatasets(ensembl_plants)
ensembl = useDataset("athaliana_eg_gene", ensembl_plants)

# Figure out what filters and attributes will be used
#filters = listFilters(ensembl)
#filters[25,]
#attributes = listAttributes(ensembl)
#attributes[5,]

# Use getBM to get attributes for values
results1 <- getBM(attributes=c('ensembl_gene_id', 'description','tair_symbol', 'name_1006', 'definition_1006'),
                  filters = 'tair_locus_model',           
                  values = atha, 
                  mart = ensembl)

# Combine with vri genes
DOM.710.ann <- cbind(matchF,results1[match(matchF$Athaliana.ortholog,results1$ensembl_gene_id),])
rm(datasets, ensembl, ensembl_plants, results1, atha)

### DOM.711
# Get A. thaliana orthologs for differentially expressed genes in domatia between genotypes
match <- merge(vriatha1,DOM.711, by=0, all="T")
match2 <- match[!is.na(match$baseMean),]
row.names(match2) <- match2$Row.names
match3 <- match2[c(2, 3, 14, 15, 16, 17, 18, 19)]
match3$Gene <- row.names(match3)
matchF <- match3[c(9, 2, 3, 4, 5, 6, 7, 8)]
names(matchF)[names(matchF) == 'V2'] <- 'Athaliana.ortholog'

# Modify file for biomaRt
col2 <- matchF$Athaliana.ortholog
test <- sapply(col2,function(ii) gsub("\\..*","",ii))
matchF$Athaliana.ortholog <- test
rm(col2,test)

library(biomaRt)

# Read in genes
atha <- matchF$Athaliana.ortholog
atha <- paste(atha, ".1", sep="")

# Load Ensembl Plants into biomaRt
ensembl_plants = useMart(biomart="plants_mart", host = "plants.ensembl.org")

# Select dataset to use
datasets <- listDatasets(ensembl_plants)
ensembl = useDataset("athaliana_eg_gene", ensembl_plants)

# Figure out what filters and attributes will be used
#filters = listFilters(ensembl)
#filters[25,]
#attributes = listAttributes(ensembl)
#attributes[5,]

# Use getBM to get attributes for values
results1 <- getBM(attributes=c('ensembl_gene_id', 'description','tair_symbol', 'name_1006', 'definition_1006'),
                  filters = 'tair_locus_model',           
                  values = atha, 
                  mart = ensembl)

# Combine with vri genes
DOM.711.ann <- cbind(matchF,results1[match(matchF$Athaliana.ortholog,results1$ensembl_gene_id),])
rm(datasets, ensembl, ensembl_plants, results1, atha)

### GENO ALONE
# Get A. thaliana orthologs for differentially expressed genes in domatia between genotypes
match <- merge(vriatha1,geno2, by=0, all="T")
match2 <- match[!is.na(match$baseMean),]
row.names(match2) <- match2$Row.names
match3 <- match2[c(2, 3, 14, 15, 16, 17, 18, 19)]
match3$Gene <- row.names(match3)
matchF <- match3[c(9, 2, 3, 4, 5, 6, 7, 8)]
names(matchF)[names(matchF) == 'V2'] <- 'Athaliana.ortholog'

# Modify file for biomaRt
col2 <- matchF$Athaliana.ortholog
test <- sapply(col2,function(ii) gsub("\\..*","",ii))
matchF$Athaliana.ortholog <- test
rm(col2,test)

library(biomaRt)

# Read in genes
atha <- matchF$Athaliana.ortholog
atha <- paste(atha, ".1", sep="")

# Load Ensembl Plants into biomaRt
ensembl_plants = useMart(biomart="plants_mart", host = "plants.ensembl.org")

# Select dataset to use
datasets <- listDatasets(ensembl_plants)
ensembl = useDataset("athaliana_eg_gene", ensembl_plants)

# Figure out what filters and attributes will be used
#filters = listFilters(ensembl)
#filters[25,]
#attributes = listAttributes(ensembl)
#attributes[5,]

# Use getBM to get attributes for values
results1 <- getBM(attributes=c('ensembl_gene_id', 'description','tair_symbol', 'name_1006', 'definition_1006'),
                  filters = 'tair_locus_model',           
                  values = atha, 
                  mart = ensembl)

# Combine with vri genes
geno2.ann <- cbind(matchF,results1[match(matchF$Athaliana.ortholog,results1$ensembl_gene_id),])
rm(datasets, ensembl, ensembl_plants, results1, atha)

############################### COMBINE ANNOTATION DATA WITH RAW READS ###############################
## DOM.GENO
#Reorder raw reads and merge with annotated differentially expressed genes
temp1 <- counts[c(2,4, 6, 1, 3, 5, 8, 10, 12, 7, 9, 11)]
temp2 <-  merge(DOM.geno.ann, temp1, by='row.names', all=FALSE)
row.names(temp2) <- temp2$Row.names

#Normalize reads, reorder them, and merge with annotated differentially expressed genes
ncounts <- counts(dds, normalized=T)
colnames(ncounts) <- paste(colnames(ncounts),"normalized",sep="-")
ncounts <- as.data.frame(ncounts)
temp3 <- ncounts[c(2,4, 6, 1, 3, 5, 8, 10, 12, 7, 9, 11)]
temp4 <-  merge(temp2, temp3, by='row.names', all=FALSE)
row.names(temp4) <- temp2$Row.names

#Clean up dataframe
DOM.GENO.final <- temp4[ -c(1,2,4) ]
DOM.GENO.final <- DOM.GENO.final[order(DOM.GENO.final$padj),] #Sort so that the genes with the lowest padj are first
rm(temp1, temp2, temp3, temp4)

## DOM.710
#Reorder raw reads and merge with annotated differentially expressed genes
temp1 <- counts[c(2,4, 6, 1, 3, 5, 8, 10, 12, 7, 9, 11)]
temp2 <-  merge(DOM.710.ann, temp1, by='row.names', all=FALSE)
row.names(temp2) <- temp2$Row.names

#Normalize reads, reorder them, and merge with annotated differentially expressed genes
ncounts <- counts(dds, normalized=T)
colnames(ncounts) <- paste(colnames(ncounts),"normalized",sep="-")
ncounts <- as.data.frame(ncounts)
temp3 <- ncounts[c(2,4, 6, 1, 3, 5, 8, 10, 12, 7, 9, 11)]
temp4 <-  merge(temp2, temp3, by='row.names', all=FALSE)
row.names(temp4) <- temp2$Row.names

#Clean up dataframe
DOM.710.final <- temp4[ -c(1,2,4) ]
DOM.710.final <- DOM.710.final[order(DOM.710.final$padj),] #Sort so that the genes with the lowest padj are first
rm(temp1, temp2, temp3, temp4)

## DOM.711
#Reorder raw reads and merge with annotated differentially expressed genes
temp1 <- counts[c(2,4, 6, 1, 3, 5, 8, 10, 12, 7, 9, 11)]
temp2 <-  merge(DOM.711.ann, temp1, by='row.names', all=FALSE)
row.names(temp2) <- temp2$Row.names

#Normalize reads, reorder them, and merge with annotated differentially expressed genes
ncounts <- counts(dds, normalized=T)
colnames(ncounts) <- paste(colnames(ncounts),"normalized",sep="-")
ncounts <- as.data.frame(ncounts)
temp3 <- ncounts[c(2,4, 6, 1, 3, 5, 8, 10, 12, 7, 9, 11)]
temp4 <-  merge(temp2, temp3, by='row.names', all=FALSE)
row.names(temp4) <- temp2$Row.names

#Clean up dataframe
DOM.711.final <- temp4[ -c(1,2,4) ]
DOM.711.final <- DOM.711.final[order(DOM.711.final$padj),] #Sort so that the genes with the lowest padj are first
rm(temp1, temp2, temp3, temp4)

## GENO ALONE
#Reorder raw reads and merge with annotated differentially expressed genes
temp1 <- counts[c(2,4, 6, 1, 3, 5, 8, 10, 12, 7, 9, 11)]
temp2 <-  merge(geno2.ann, temp1, by='row.names', all=FALSE)
row.names(temp2) <- temp2$Row.names

#Normalize reads, reorder them, and merge with annotated differentially expressed genes
ncounts <- counts(dds, normalized=T)
colnames(ncounts) <- paste(colnames(ncounts),"normalized",sep="-")
ncounts <- as.data.frame(ncounts)
temp3 <- ncounts[c(2,4, 6, 1, 3, 5, 8, 10, 12, 7, 9, 11)]
temp4 <-  merge(temp2, temp3, by='row.names', all=FALSE)
row.names(temp4) <- temp2$Row.names

#Clean up dataframe
geno2.final <- temp4[ -c(1,2,4) ]
geno2.final <- geno2.final[order(geno2.final$padj),] #Sort so that the genes with the lowest padj are first
rm(temp1, temp2, temp3, temp4)

############################### SAVE OUTPUT FILE ###############################
#DOM.GENO
#write.csv(DOM.GENO.final, file="DE_genes_Domatia_710_v_711.csv", row.names=FALSE)

# DOM.710
#write.csv(DOM.710.final, file="DE_genes_Domatia_V_Leaf_588710.csv", row.names=FALSE)

# DOM.711
#write.csv(DOM.711.final, file="DE_genes_Domatia_V_Leaf_588711.csv", row.names=FALSE)

# GENO ONLY
#write.csv(geno2.final, file="DE_genes_Control_710_v_711.csv", row.names=FALSE)

############################### WHICH GENES DIFFER BETWEEN LISTS? ###############################

temp1 <- DOM.710.final[!rownames(DOM.710.final)%in%rownames(DOM.711.final),]

temp2 <- DOM.711.final[!rownames(DOM.711.final)%in%rownames(DOM.710.final),]


############################### WHICH GENES ARE SHARED BETWEEN LISTS? ###############################

# What genes unique to the domatia 710 and 711 lists are differentially expressed
# between genotypes?

temp5 <- temp1[rownames(temp1)%in%rownames(geno2.final),]

temp6 <- temp2[rownames(temp2)%in%rownames(geno2.final),]

# What genes in the domatia 710 and 711 lists are differentially expressed
# in both genotypes?

dom.both <- DOM.710.final[rownames(DOM.710.final)%in%rownames(DOM.711.final),]

# DIFFERENT VERSION BELOW, READING IN FILES ALREADY MADE
concon <- read.csv("differentially_expressed_genes/DE_genes_Control_710_v_711.csv")
row.names(concon) <- concon$Gene

domdom <- read.csv("differentially_expressed_genes/DE_genes_Domatia_710_v_711.csv")
row.names(domdom) <- domdom$Gene

domcon.710 <- read.csv("differentially_expressed_genes/DE_genes_Domatia_V_Leaf_588710.csv")
row.names(domcon.710) <- domcon.710$Gene

domcon.711 <- read.csv("differentially_expressed_genes/DE_genes_Domatia_V_Leaf_588711.csv")
row.names(domcon.711) <- domcon.711$Gene

# Shared between domatia and domatia
a <- domcon.710[rownames(domcon.710)%in%rownames(domcon.711),]

# Shared between control and domatia 710
b <- concon[rownames(concon)%in%rownames(domcon.710),]

# Shared between control and domatia 711
c <- concon[rownames(concon)%in%rownames(domcon.711),]

# Shared between domatia and domatia 710
d <- domcon.710[rownames(domcon.710)%in%rownames(domdom),]

# Shared between domatia and domatia 711
e <- domcon.711[rownames(domcon.711)%in%rownames(domdom),]

# Shared between control and domatia
f <- concon[rownames(concon)%in%rownames(domdom),]