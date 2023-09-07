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
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)

######################## LOAD DATA ########################
# LOAD DIFFERENTIALLY EXPRESSED GENES
concon <- read.csv("differentially_expressed_genes/DE_genes_Control_710_v_711.csv")
row.names(concon) <- concon$Gene
domdom <- read.csv("differentially_expressed_genes/DE_genes_Domatia_710_v_711.csv")
row.names(domdom) <- domdom$Gene
domcon.710 <- read.csv("differentially_expressed_genes/DE_genes_Domatia_V_Leaf_588710.csv")
row.names(domcon.710) <- domcon.710$Gene
domcon.711 <- read.csv("differentially_expressed_genes/DE_genes_Domatia_V_Leaf_588711.csv")
row.names(domcon.711) <- domcon.711$Gene

# LOAD ARABIDOPSIS DATA AND GENES IF GETTING GO TERMS
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


######################## COMPARING COMPLEX HEATMAP TO OTHER OPTIONS ########################
# Make test matrix
genes.df <- domdom[c(23,24,25,26,27,28,29,30,31,32,33,34)]
colnames(genes.df) <- c("588710 1 Control", "588710 2 Control", "588710 3 Control", 
                        "588710 1 Domatia", "588710 2 Domatia", "588710 3 Domatia" ,
                        "588711 1 Control", "588711 2 Control", "588711 3 Control",
                        "588711 1 Domatia", "588711 2 Domatia", "588711 3 Domatia")
genes.mat <- as.matrix(genes.df)

# Heatmap with image function
image(apply(genes.mat,1,scale))


# Making sure that I am normalizing how base R's heatmap normalizes data
tmp<-t(apply(genes.mat,1,scale))
col_fun<-colorRamp2(seq(min(tmp),max(tmp),length.out=256),
                    colorRampPalette(c("blue","white","red"))(256))
Heatmap(tmp, cluster_rows = FALSE, cluster_columns = FALSE,col=col_fun)

col<- colorRampPalette(c("blue", "white", "red"))(256)
heatmap(genes.mat, col=col, Rowv=NA, Colv=NA)

# Complex heatmap seems to be the best option! I can normalize it like base R's heatmap,
# but make it look nice.

######################## COMPARING LOG TRANSFORMED VERSUS NON-TRANSFORMED ########################
# Plotting log transformed and non-transformed in case Chad or Marge have an issue
# with the transformed version.

# Make test matrix
genes.df <- domdom[c(23,24,25,26,27,28,29,30,31,32,33,34)]
colnames(genes.df) <- c("588710 1 Control", "588710 2 Control", "588710 3 Control", 
                        "588710 1 Domatia", "588710 2 Domatia", "588710 3 Domatia" ,
                        "588711 1 Control", "588711 2 Control", "588711 3 Control",
                        "588711 1 Domatia", "588711 2 Domatia", "588711 3 Domatia")
genes.mat <- as.matrix(genes.df)

# LOG TRANSFORMED
logt <- log(genes.mat+0.01)
temp1 <- t(apply(logt, 1, scale))
col_fun<-colorRamp2(seq(min(temp1),max(temp1),length.out=256),
                    colorRampPalette(c("blue","white","red"))(256))
Heatmap(temp1, cluster_rows = FALSE, cluster_columns = FALSE,col=col_fun)

# WITHOUT LOG TRANSFORMATION
temp2 <- t(apply(genes.mat, 1, scale))
col_fun<-colorRamp2(seq(min(temp2),max(temp2),length.out=256),
                    colorRampPalette(c("blue","white","red"))(256))
Heatmap(temp2, cluster_rows = FALSE, cluster_columns = FALSE,col=col_fun)

######################## MAKE FIGURE 3 ########################
cellwall <- read.csv("GO-term-enrichment/DE_genes_Control_V_Domatia_588710_cell_wall_gene_list.csv",
                     sep="\t", row.names=NULL, header = FALSE)
genes.cellwall <- cellwall$V2
genes.df <- domcon.710[domcon.710$ensembl_gene_id %in% genes.cellwall,]
genes.df1 <- genes.df[genes.df$Gene %in% domcon.711$Gene,] # Run if we only want genes shared between genotypes
genes.df2 <- genes.df1[c(23,24,25,26,27,28,29,30,31,32,33,34)] # Modify if code above is run
colnames(genes.df2) <- c("588710 1 Control", "588710 2 Control", "588710 3 Control", 
                        "588710 1 Domatia", "588710 2 Domatia", "588710 3 Domatia" ,
                        "588711 1 Control", "588711 2 Control", "588711 3 Control",
                        "588711 1 Domatia", "588711 2 Domatia", "588711 3 Domatia")
genes.mat <- as.matrix(genes.df2)
logt <- log(genes.mat+0.01)
temp1 <- t(apply(logt, 1, scale))
col_fun<-colorRamp2(seq(min(temp1),max(temp1),length.out=256),
                    colorRampPalette(c("blue","white","red"))(256))
Heatmap(temp1, cluster_rows = FALSE, cluster_columns = FALSE,col=col_fun)

######################## MAKE FIGURE 4 ########################

genes.meta <- c("",
                )

######################## MAKE FIGURE 5 ########################
auxin <- read.csv("GO-term-enrichment/for-heatmaps/auxin/588710-auxin.txt",
                     sep="\t", row.names=NULL, header = FALSE)
auxin.genes <- c("AT5G57390", #AIL5 #Start of auxin synthesis regulators 
                 "AT5G10510", #AIL6
                 "AT5G54510", #GH3.6 #Start of auxin synthases
                 "AT1G28130", #GH3.17
                 "AT1G73590", #PIN1 #Start of auxin transporters
                 "AT1G70940", #PIN3
                 "AT1G75500", #WAT1
                 "AT3G23050", #IAA7 #Start of transcriptional regulation via auxin signaling
                 "AT2G46990", #IAA20
                 "AT2G33310", #IAA13
                 "AT3G16500", #IAA26
                 "AT5G43700", #IAA4
                 "AT3G15540", #IAA19
                 "AT5G57420", #IAA33
                 "AT1G04550", #IAA12
                 "AT5G65670", #IAA9
                 "AT4G29080", #IAA27
                 "AT1G19850", #ARF5
                 "AT4G34810", #SAUR5 #Start of genes upregulated by auxin
                 "AT4G34800", #SAUR4
                 "AT2G21210", #SAUR6
                 "AT3G12955", #SAUR74
                 "AT1G75590", #SAUR52
                 "AT4G38840", #SAUR14
                 "AT1G29450", #SAUR64
                 "AT1G29510", #SAUR67
                 "AT1G29420" #SAUR61
                 )

genes.df <- domcon.710[domcon.710$ensembl_gene_id %in% auxin.genes,]
genes.df1 <- genes.df[genes.df$Gene %in% domcon.711$Gene,] # Run if we only want genes shared between genotypes
genes.df1 <- genes.df1[order(match(genes.df1$ensembl_gene_id,auxin.genes)),]
genes.df2 <- genes.df1[c(23,24,25,26,27,28,29,30,31,32,33,34)] # Modify if code above is run
colnames(genes.df2) <- c("588710 1 Control", "588710 2 Control", "588710 3 Control", 
                         "588710 1 Domatia", "588710 2 Domatia", "588710 3 Domatia" ,
                         "588711 1 Control", "588711 2 Control", "588711 3 Control",
                         "588711 1 Domatia", "588711 2 Domatia", "588711 3 Domatia")

rownames(genes.df2) <- c("GH3.6", "GH3.17-like", "Auxin efflux carrier component 1c",
  "Auxin efflux carrier component 3-like", "WAT1", "IAA14", "IAA13-like", "Auxin-induced protein 22D-like (1)",
  "Auxin-induced protein 22D-like (2)", "AUX22-like", "IAA9-like", "IAA27",
  "LOC117913455", "SAUR21-like (1)", "SAUR21-like (2)", "SAUR21-like (3)", "X15-like",
  "SAUR64-like", "SAUR67-like", "SAUR68-like")

genes.mat <- as.matrix(genes.df2)

# #Log transformed values below - I think it looks better not transformed
# logt <- log(genes.mat+0.001)
# temp1 <- t(apply(logt, 1, scale))
# col_fun<-colorRamp2(seq(min(temp1),max(temp1),length.out=256),
#                     colorRampPalette(c("blue","white","red"))(256))
# Heatmap(temp1, cluster_rows = FALSE, cluster_columns = FALSE,col=col_fun)

temp1 <- t(apply(genes.df2, 1, scale))
col_fun<-colorRamp2(seq(min(temp1),max(temp1),length.out=256),
                    colorRampPalette(c("blue","white","red"))(256))
Heatmap(temp1, cluster_rows = FALSE, cluster_columns = FALSE,col=col_fun)

column_split = rep("Control 588710", 12)
column_split[4:6] = "Domatia 588710"
column_split[7:9] = "Control 588711"
column_split[10:12] = "Domatia 588711"

row_split = rep("Auxin synthases", 20)
row_split[3:5] = "Auxin transporters"
row_split[6:12] = "Transcriptional regulation\n via auxin signaling"
row_split[13:20] = "Genes upregulated by auxin"

#Legend to right
Heatmap(temp1, 
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        col=col_fun,
        heatmap_legend_param = list(
          title = "Z-Score", 
          border = "black", 
          legend_height = unit(6, "cm"),
          title_gp = gpar(size = 14, fontface = "bold")),
        column_split = factor(column_split, levels = c("Control 588710", "Domatia 588710", "Control 588711", "Domatia 588711")),
        cluster_column_slices = FALSE,
        column_gap = unit(2, "mm"),
        border = TRUE,
        row_split = factor(row_split, levels = c("Auxin synthases", "Auxin transporters", "Transcriptional regulation\n via auxin signaling", "Genes upregulated by auxin")),
        row_title_rot = 0,
        row_gap = unit(2, "mm")
)
#Legend on bottom
htmp = Heatmap(temp1, 
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        col=col_fun,
        heatmap_legend_param = list(
          title = "Z-Score", 
          border = "black", 
          legend_width = unit(6, "cm"),
          title_gp = gpar(size = 14, fontface = "bold"),
          direction = "horizontal"),
        column_split = factor(column_split, levels = c("Control 588710", "Domatia 588710", "Control 588711", "Domatia 588711")),
        cluster_column_slices = FALSE,
        column_gap = unit(2, "mm"),
        border = TRUE,
        row_split = factor(row_split, levels = c("Auxin synthases", "Auxin transporters", "Transcriptional regulation\n via auxin signaling", "Genes upregulated by auxin")),
        row_title_rot = 0,
        row_gap = unit(2, "mm"),
        row_names_gp = gpar(fontsize = 11)
)

draw(htmp, heatmap_legend_side="bottom")
######################## MAKE FIGURE 8 ########################
dom.genes <- c("LOC117934313",
               "LOC117912178",
               "LOC117912293",
               "LOC117933133",
               "LOC117930032",
               "LOC117906993",
               "LOC117927833",
               "LOC117904283",
               "LOC117922028") #Still need to add random genes

genes.df <- DOM.GENO.final[DOM.GENO.final$Gene %in% dom.genes,]
#genes.df1 <- genes.df[genes.df$Gene %in% domcon.711$Gene,] # Run if we only want genes shared between genotypes
genes.df1 <- genes.df[order(match(genes.df$Gene,dom.genes)),]
genes.df2 <- genes.df1[c(25,26,27,28,29,30,31,32,33,34, 35, 36)] # Modify if code above is run
colnames(genes.df2) <- c("588710 1 Control", "588710 2 Control", "588710 3 Control", 
                         "588710 1 Domatia", "588710 2 Domatia", "588710 3 Domatia" ,
                         "588711 1 Control", "588711 2 Control", "588711 3 Control",
                         "588711 1 Domatia", "588711 2 Domatia", "588711 3 Domatia")

genes.mat <- as.matrix(genes.df2)
logt <- log(genes.mat+0.01)
temp1 <- t(apply(logt, 1, scale))
col_fun<-colorRamp2(seq(min(temp1),max(temp1),length.out=256),
                    colorRampPalette(c("blue","white","red"))(256))
Heatmap(temp1, cluster_rows = FALSE, cluster_columns = FALSE,col=col_fun)

####################### GETTING VRIPARIA PRODUCT NAMES ########################
gff.id <- read.csv("vriparia-genomic-id.csv", sep=";", header = FALSE)

test <- as.matrix(gff.id)
genes<-matrix(grepl("^gene=",test),nrow(gff.id),ncol(gff.id))
any(rowSums(genes)>1) #gene specified only once per row
tmp<-which(genes,arr.ind=TRUE)
genes<-gsub("^gene=","",test[genes])
genes<-genes[match(1:nrow(test),tmp[,1])]
products<-matrix(grepl("^product=",test),nrow(gff.id),ncol(gff.id))
any(rowSums(products)>1) #products specified only once per row
tmp<-which(products,arr.ind=TRUE)
products<-gsub("^product=","",test[products])
products<-products[match(1:nrow(test),tmp[,1])]
dict <- as.data.frame(unique(cbind(genes,products)))

temp1 <- na.omit(dict[dict$genes %in% dom.genes,])
temp2 <- temp1[!duplicated(temp1$genes),]

temp1 <- na.omit(dict[dict$genes %in% rownames(genes.df2),])
temp2 <- temp1[!duplicated(temp1$genes),]
