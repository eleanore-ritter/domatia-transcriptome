# Set working directory and load necessary packages
setwd("C:/Users/rittere5/OneDrive - Michigan State University/Vitis-domatia/")
#setwd("C:/Users/elean/OneDrive - Michigan State University/Vitis-domatia/")
library(DESeq2)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(ggh4x)
library(reshape2)
library(gplots)

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
dds = DESeqDataSetFromMatrix(countData = counts, colData = colData,  design = ~treatment)
#Make sure Control is the first level of treatment
dds$treatment <- relevel(dds$treatment, ref = "Control")
#Make sure that dds is properly formatted
as.data.frame(colData(dds))

##Run DESeq
dds <- DESeq(dds)
res <- results(dds)

############################### EXTRACTING AND FILTERING SPECIFIC SAMPLES ################################
## Double check result names
resultsNames(dds) #Should only have "intercept" and "treatment_X_vs_Y"

## Filter differentially expressed genes
deseq1 <- results(dds)
df.deseq1 <- as.data.frame(deseq1)
res1.deseq1<-df.deseq1[(df.deseq1$log2FoldChange>1|df.deseq1$log2FoldChange< -1)& !is.na(df.deseq1$log2FoldChange),] #Only keep genes with a log2FoldChange of at least +/-1
final.deseq <- res1.deseq1[(res1.deseq1$padj<0.05)&!is.na(res1.deseq1$padj),] #Only keep genes with a padj below 0.05
final.deseq <- final.deseq[order(final.deseq$padj),] #Sort so that the genes with the lowest padj are first

############################### ATHA ORTHOLOGS ###############################
# Load diamond tsv file for Vitis riparia and Arabidopsis thaliana
vriatha <- read.csv("Vriparia-Athaliana.tsv", header = FALSE, sep = "\t")
vriatha1 <- vriatha %>% distinct(V1, .keep_all = TRUE)
row.names(vriatha1) <- vriatha1$V1

# Get A. thaliana orthologs for differentially expressed genes
match <- merge(vriatha1,final.deseq, by=0, all="T")
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
final.deseq.ann <- cbind(matchF,results1[match(matchF$Athaliana.ortholog,results1$ensembl_gene_id),])
final.deseq.ann <- final.deseq.ann[order(final.deseq.ann$padj),] #Sort so that the genes with the lowest padj are first

## Add in GO or PO terms if wanted
# results2 <- getBM(attributes=c('ensembl_gene_id', 'description','tair_symbol', 'go_id','name_1006', 'definition_1006', 'po_id', 'po_name_1006'),
#                   filters = 'tair_locus_model',           
#                   values = atha, 
#                   mart = ensembl)
# final.deseq.ann2 <- cbind(matchF,results2[match(matchF$Athaliana.ortholog,results2$ensembl_gene_id),])
# final.deseq.ann2 <- final.deseq.ann2[order(final.deseq.ann2$padj),]

rm(datasets, ensembl, ensembl_plants, results1, atha)

############################### COMBINE ANNOTATION DATA WITH READS ###############################
#Reorder raw reads and merge with annotated differentially expressed genes
temp1 <- counts[c(2,4, 6, 1, 3, 5, 8, 10, 12, 7, 9, 11)]
temp2 <-  merge(final.deseq.ann, temp1, by='row.names', all=FALSE)
row.names(temp2) <- temp2$Row.names

#Normalize reads, reorder them, and merge with annotated differentially expressed genes
ncounts <- counts(dds, normalized=T)
colnames(ncounts) <- paste(colnames(ncounts),"normalized",sep="-")
ncounts <- as.data.frame(ncounts)
temp3 <- ncounts[c(2,4, 6, 1, 3, 5, 8, 10, 12, 7, 9, 11)]
temp4 <-  merge(temp2, temp3, by='row.names', all=FALSE)
row.names(temp4) <- temp2$Row.names

#Clean up dataframe
final <- temp4[ -c(1,2,4) ]
final <- final[order(final$padj),] #Sort so that the genes with the lowest padj are first

############################### SAVE OUTPUT FILE ###############################

#write.csv(final, file="DE_genes_Domatia_v_Leaf.csv", row.names = FALSE)

############################### HEATMAP ###############################
final.deseq <- final.deseq[order(final.deseq$padj),]
genes <- row.names(final.deseq)
genes.df <- counts[row.names(counts) %in% genes,] #Could use normalized counts or raw - normalized is probably better, but I need to look into it
genes.df <- genes.df[c(2,4, 6, 1, 3, 5, 8, 10, 12, 7, 9, 11)]
genes.mat <- as.matrix(genes.df)
heatmap(genes.mat, Rowv = NA, Colv = NA)

#Plot only a subset
final.deseq <- final.deseq[order(final.deseq$padj),]
final.deseq.sub <- head(final.deseq, 50)
genes <- row.names(final.deseq.sub)
genes.df <- ncounts[row.names(ncounts) %in% genes,] #Could use normalized counts or raw - normalized is probably better, but I need to look into it
genes.df <- genes.df[c(2,4, 6, 1, 3, 5, 8, 10, 12, 7, 9, 11)]
genes.mat <- as.matrix(genes.df)
heatmap(genes.mat, Rowv = NA, Colv = NA)

############################### MAKE GENE EXPRESSION PLOTS ###############################
# Load data
#final <- read.csv("DE_genes_Domatia_v_Leaf.csv")
#row.names(final) <- final$Gene

# Restructure data
genename <- c("LOC117907723") # Define gene of interest
gene <- final[rownames(final)%in%genename,23:34] # Get normalized reads for only the gene of interest
gene.long <- pivot_longer(gene, everything()) # Transpose from wide to long format for plotting
gene.long$sample <- c("1", "4", "5", "1", "4", "5", "6", "8", "9", "6", "8", "9")
gene.long$treatment <- c("Control", "Control", "Control", "Domatia", "Domatia", "Domatia",
                         "Control", "Control", "Control", "Domatia", "Domatia", "Domatia")
gene.long$genotype <- c("588710", "588710","588710", "588710", "588710","588710",
                        "588711","588711","588711","588711","588711","588711")

# Plotting
# domcolor <- c("#114b5f")
# concolor <- c("#1a936f")

domcolor <- c("#558B6E")
concolor <- c("#383B56")


ggplot(gene.long, aes(x=sample, y=value)) +
  geom_point(aes(color=treatment), size=3) +
  scale_color_manual(values = c(concolor, domcolor)) +
  theme_classic() +
  labs(x = "", y = "Normalized Read Counts") +
  scale_x_discrete(labels=c("1", "2", "3", "1", "2", "3")) +
  annotate(geom = "text", x = 2 + 3 * (0:1), y = -2500, label = unique(gene.long$genotype), size = 4) +
  annotate("segment", x = 0.75, xend = 3.25, y = -2050, yend = -2050) +
  annotate("segment", x = 3.75, xend = 6.25, y = -2050, yend = -2050) +
  coord_cartesian(ylim = c(0, 16000), clip = "off") +
  theme(axis.title.x=element_text(colour="black", size=12),
        axis.title.y=element_text(colour="black", size=12),
        axis.text.x=element_text(colour="black", size=10),
        axis.text.y = element_text(colour="black", size=10),
        legend.position = "right",
        legend.justification = "center",
        legend.title = element_blank(),
        legend.key.size = unit(0.25, "cm"),
        legend.text = element_text(colour="black", size=10))

############################### HEATMAPS FOR DIFFERENT GENE GROUPS ###############################
#Plot only a subset
cell.wall <- remove_missing(final[final$name_1006=="cell wall organization",])
genes <- rownames(cell.wall)
genes.df <- ncounts[row.names(ncounts) %in% genes,] #Could use normalized counts or raw - normalized is probably better, but I need to look into it
genes.df <- genes.df[c(2,4, 6, 1, 3, 5, 8, 10, 12, 7, 9, 11)]
colnames(genes.df) <- c("588710 1 Control", "588710 2 Control", "588710 3 Control", 
                        "588710 1 Domatia", "588710 2 Domatia", "588710 3 Domatia" ,
                        "588711 1 Control", "588711 2 Control", "588711 3 Control",
                        "588711 1 Domatia", "588711 2 Domatia", "588711 3 Domatia")
genes.mat <- as.matrix(genes.df)
heatmap(genes.mat, Rowv = NA, Colv = NA)

#Plot with ggplot2
heatmap.2(genes.mat, scale="row", dendrogram = c("none"))
