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

overlap <- domcon.710[domcon.710$Gene %in% domcon.711$Gene,]
writeClipboard(paste0(as.vector(na.omit(gsub("\\.[0-9]", "", overlap$ensembl_gene_id)))))
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

cellwall.genes <- c("AT5G13870", #XTH5 #Start of xyloglucan synthesis and modification
                    "AT2G14620", #XTH10
                    "AT4G14130", #XTH15
                    "AT4G25810", #XTH23
                    "AT2G36870", #XTH32
                    "AT5G20950", #BGLC1
                    "AT1G19300", #GATL1 #Start of xylan synthesis and modification
                    "AT3G50220", #IRX15 
                    "AT3G18660", #GUX1
                    "AT5G54690", #GAUT12
                    "AT1G09610", #GXM1
                    "AT1G33800", #GXM3
                    "AT1G73140", #TBL31
                    "AT2G38320", #TBL34
                    "AT3G05620", #PME22 #Start of pectin biosynthesis and modification
                    "AT3G43270", #PME32
                    "AT1G23200", #PME6
                    "AT2G26440", #PME12
                    "AT1G11580", #PME18
                    "AT5G19730", #PME53
                    "AT1G05310", #PME8
                    "AT5G63180", #PLL15
                    "AT1G67750", #PLL16
                    "AT3G07010", #PLL20
                    "AT3G29360", #UGD2
                    "AT5G03260", #LAC11 #Start of lignin biosynthesis and modification
                    "AT5G60020", #LAC17
                    "AT2G38080", #LAC4
                    "AT4G36220", #FAH1
                    "AT4G34050", #CCOAOMT1
                    "AT1G67980", #CCOAMT
                    "AT5G05340", #PER52
                    "AT4G01070" #UGT72B1
                    )

genes.df <- domcon.710[domcon.710$ensembl_gene_id %in% cellwall.genes,]
genes.df1 <- genes.df[genes.df$Gene %in% domcon.711$Gene,] # Run if we only want genes shared between genotypes
genes.df1 <- genes.df1[order(match(genes.df1$ensembl_gene_id,cellwall.genes)),]
genes.df2 <- genes.df1[c(23,24,25,26,27,28,29,30,31,32,33,34)] # Modify if code above is run
colnames(genes.df2) <- c("588710 1 Control", "588710 2 Control", "588710 3 Control", 
                         "588710 1 Domatia", "588710 2 Domatia", "588710 3 Domatia" ,
                         "588711 1 Control", "588711 2 Control", "588711 3 Control",
                         "588711 1 Domatia", "588711 2 Domatia", "588711 3 Domatia")

rownames(genes.df2) <- c("Probable xyloglucan endotransglucosylase/\nhydrolase protein B",
                         "Probable xyloglucan endotransglucosylase/\nhydrolase protein 10",
                         "Xyloglucan endotransglucosylase/\nhydrolase 2-like",
                         "Pectinesterase",
                         "Glucuronoxylan\n4-O-methyltransferase 3",
                         "Probable pectinesterase 53",
                         "Pectate lyase-like",
                         "Probable pectate lyase 5",
                         "Laccase-11-like",
                         "Laccase-17-like (1)",
                         "Laccase-17-like (2)",
                         "Laccase-17-like (3)",
                         "Laccase-17-like (4)",
                         "Laccase-17-like (5)",
                         "Hydroquinone\nglucosyltransferase-like")

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

column_split = rep("Control\n 588710", 12)
column_split[4:6] = "Domatia\n 588710"
column_split[7:9] = "Control\n 588711"
column_split[10:12] = "Domatia\n 588711"

row_split = rep("Xyloglucan\nbiosynthesis\nand modification", 15)
row_split[4] = "Xylan biosynthesis\n and modification"
row_split[6:8] = "Pectin biosynthesis\n and modification"
row_split[9:15] = "Lignin biosynthesis\n and modification"

#Legend to right
# Heatmap(temp1, 
#         cluster_rows = FALSE, 
#         cluster_columns = FALSE,
#         col=col_fun,
#         heatmap_legend_param = list(
#           title = "Z-Score", 
#           border = "black", 
#           legend_height = unit(6, "cm"),
#           title_gp = gpar(size = 14, fontface = "bold")),
#         column_split = factor(column_split, levels = c("Control 588710", "Domatia 588710", "Control 588711", "Domatia 588711")),
#         cluster_column_slices = FALSE,
#         column_gap = unit(2, "mm"),
#         border = TRUE,
#         row_split = factor(row_split, levels = c("Auxin synthases", "Auxin transporters", "Transcriptional regulation\n via auxin signaling", "Genes upregulated\n by auxin")),
#         row_title_rot = 0,
#         row_gap = unit(2, "mm")
# )

#Legend on bottom
htmp = Heatmap(temp1, 
               cluster_rows = FALSE, 
               cluster_columns = FALSE,
               col=col_fun,
               heatmap_legend_param = list(
                 title = "Z-Score", 
                 border = "black", 
                 legend_width = unit(6, "cm"),
                 title_gp = gpar(size = 12, fontface = "bold"),
                 direction = "horizontal"),
               column_split = factor(column_split, levels = c("Control\n 588710", "Domatia\n 588710", "Control\n 588711", "Domatia\n 588711")),
               cluster_column_slices = FALSE,
               column_gap = unit(2, "mm"),
               border = TRUE,
               row_split = factor(row_split, levels = c("Xyloglucan\nbiosynthesis\nand modification", "Xylan biosynthesis\n and modification", 
                                                        "Pectin biosynthesis\n and modification", "Lignin biosynthesis\n and modification")),
               row_title_rot = 0,
               row_gap = unit(2, "mm"),
               row_names_max_width = max_text_width(
                 rownames(temp1), 
                 gp = gpar(fontsize = 10)
               ),
               row_names_gp = gpar(fontsize = 10),
               row_title_gp = gpar(fontsize = 12, fontface = "bold")
)


draw(htmp, heatmap_legend_side="bottom")

######################## MAKE FIGURE 4 ########################

meta.genes <- c("AT5G09220", #AAP2 #Start of amino acid transport
                "AT1G77380", #AAP3
                "AT1G10010", #AAP8
                "AT5G04770", #CAT6
                "AT2G39510", #UMAMIT14
                "AT3G28960", #hexon #Start of carb metabolism
                "AT1G73370", #SUS6
                "AT1G12240", #VIN2
                "AT1G47840", #HXK3
                "AT4G10260", #FRK4
                "AT5G49690", #UGT91C1
                "AT1G22360", #AtUGT85A2
                "AT1G22400", #UGT85A1
                "AT4G01070", #GT72B1
                "AT2G43820", #GT
                "AT5G40390", #RFS5 
                "AT4G15920", #SWEET17 #Start of carb transport
                "AT1G22710", #SUC2
                "AT1G35720" #ANN1
)

genes.df <- domcon.710[domcon.710$ensembl_gene_id %in% meta.genes,]
genes.df1 <- genes.df[genes.df$Gene %in% domcon.711$Gene,] # Run if we only want genes shared between genotypes
genes.df1 <- genes.df1[order(match(genes.df1$ensembl_gene_id,meta.genes)),]
genes.df2 <- genes.df1[c(23,24,25,26,27,28,29,30,31,32,33,34)] # Modify if code above is run
colnames(genes.df2) <- c("588710 1 Control", "588710 2 Control", "588710 3 Control", 
                         "588710 1 Domatia", "588710 2 Domatia", "588710 3 Domatia" ,
                         "588711 1 Control", "588711 2 Control", "588711 3 Control",
                         "588711 1 Domatia", "588711 2 Domatia", "588711 3 Domatia")

rownames(genes.df2) <- c("Amino acid permease 3-like (1)",
                         "Amino acid permease 3-like (2)",
                         "Amino acid permease 8-like",
                         "Cationic amino acid transporter 6",
                         "WAT1-related protein",
                         "Amino acid transporter AVT1I-like (1)",
                         "Amino acid transporter AVT1I-like (2)",
                         "Sucrose synthase 7-like",
                         "Acid beta-fructofuranosidase",
                         "Hexokinase-2",
                         "Probable fructokinase-5",
                         "UDP-glycosyltransferase 91C1-like",
                         "7-deoxyloganetin glucosyltransferase-like",
                         "7-deoxyloganetic acid glucosyltransferase-like",
                         "Hydroquinone glucosyltransferase-like",
                         "UDP-glycosyltransferase 74F2-like",
                         "Galactinol--sucrose galactosyltransferase-like (1)",
                         "Galactinol--sucrose galactosyltransferase-like (2)",
                         "Bidirectional sugar transporter SWEET17",
                         "Sucrose transport protein SUC2-like",
                         "Annexin D1",
                         "Annexin D4-like")

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

column_split = rep("Control\n 588710", 12)
column_split[4:6] = "Domatia\n 588710"
column_split[7:9] = "Control\n 588711"
column_split[10:12] = "Domatia\n 588711"

row_split = rep("Amino acid\ntransport", 22)
row_split[8:18] = "Carbohydrate\nmetabolism"
row_split[19:22] = "Carbohydrate\ntransport"

#Legend on bottom
htmp = Heatmap(temp1, 
               cluster_rows = FALSE, 
               cluster_columns = FALSE,
               col=col_fun,
               heatmap_legend_param = list(
                 title = "Z-Score", 
                 border = "black", 
                 legend_width = unit(6, "cm"),
                 title_gp = gpar(size = 12, fontface = "bold"),
                 direction = "horizontal"),
               column_split = factor(column_split, levels = c("Control\n 588710", "Domatia\n 588710", "Control\n 588711", "Domatia\n 588711")),
               cluster_column_slices = FALSE,
               column_gap = unit(2, "mm"),
               border = TRUE,
               row_split = factor(row_split, levels = c("Amino acid\ntransport", "Carbohydrate\nmetabolism", 
                                                        "Carbohydrate\ntransport")),
               row_title_rot = 0,
               row_gap = unit(2, "mm"),
               row_names_max_width = max_text_width(
                 rownames(temp1), 
                 gp = gpar(fontsize = 10)
               ),
               row_names_gp = gpar(fontsize = 10),
               row_title_gp = gpar(fontsize = 12, fontface = "bold")
)


draw(htmp, heatmap_legend_side="bottom")

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

column_split = rep("Control\n 588710", 12)
column_split[4:6] = "Domatia\n 588710"
column_split[7:9] = "Control\n 588711"
column_split[10:12] = "Domatia\n 588711"

row_split = rep("Auxin synthases", 20)
row_split[3:5] = "Auxin transporters"
row_split[6:12] = "Transcriptional regulation\n via auxin signaling"
row_split[13:20] = "Genes upregulated\n by auxin"

#Legend to right
# Heatmap(temp1, 
#         cluster_rows = FALSE, 
#         cluster_columns = FALSE,
#         col=col_fun,
#         heatmap_legend_param = list(
#           title = "Z-Score", 
#           border = "black", 
#           legend_height = unit(6, "cm"),
#           title_gp = gpar(size = 14, fontface = "bold")),
#         column_split = factor(column_split, levels = c("Control 588710", "Domatia 588710", "Control 588711", "Domatia 588711")),
#         cluster_column_slices = FALSE,
#         column_gap = unit(2, "mm"),
#         border = TRUE,
#         row_split = factor(row_split, levels = c("Auxin synthases", "Auxin transporters", "Transcriptional regulation\n via auxin signaling", "Genes upregulated\n by auxin")),
#         row_title_rot = 0,
#         row_gap = unit(2, "mm")
# )

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
        column_split = factor(column_split, levels = c("Control\n 588710", "Domatia\n 588710", "Control\n 588711", "Domatia\n 588711")),
        cluster_column_slices = FALSE,
        column_gap = unit(2, "mm"),
        border = TRUE,
        row_split = factor(row_split, levels = c("Auxin synthases", "Auxin transporters", "Transcriptional regulation\n via auxin signaling", "Genes upregulated\n by auxin")),
        row_title_rot = 0,
        row_gap = unit(2, "mm"),
        row_names_gp = gpar(fontsize = 11)
)

draw(htmp, heatmap_legend_side="bottom")
######################## MAKE FIGURE 6 ########################

biotic.genes <- c("AT1G12220", #NBS-LRR class
                  "AT3G51570",
                  "AT5G17680",
                  "AT5G36930",
                  "AT3G61460", #BRH1 #Chitin
                  "AT5G63380", #JA and methyljasmonate biosynthesis #JA signaling
                  "AT1G20510",
                  "AT4G08850",
                  "AT1G19640",
                  "AT2G41370", #JA signaling
                  "AT1G56650",
                  
)

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

# temp1 <- na.omit(dict[dict$genes %in% dom.genes,])
# temp2 <- temp1[!duplicated(temp1$genes),]

temp1 <- na.omit(dict[dict$genes %in% rownames(genes.df2),])
temp2 <- temp1[!duplicated(temp1$genes),]
