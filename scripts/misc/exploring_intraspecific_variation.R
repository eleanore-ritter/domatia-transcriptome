# Set working directory and load necessary packages
#setwd("C:/Users/rittere5/OneDrive - Michigan State University/Vitis-domatia/")
setwd("C:/Users/elean/OneDrive - Michigan State University/Vitis-domatia/")

library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(ggcorrplot)
library(viridis)
library(topGO)

############################# LOAD AND SET UP DATA #############################

concon <- read.csv("differentially_expressed_genes/DE_genes_Control_710_v_711.csv")
row.names(concon) <- concon$Gene

domdom <- read.csv("differentially_expressed_genes/DE_genes_Domatia_710_v_711.csv")
row.names(domdom) <- domdom$Gene

domcon.710 <- read.csv("differentially_expressed_genes/DE_genes_Domatia_V_Leaf_588710.csv")
row.names(domcon.710) <- domcon.710$Gene

domcon.711 <- read.csv("differentially_expressed_genes/DE_genes_Domatia_V_Leaf_588711.csv")
row.names(domcon.711) <- domcon.711$Gene

# Get genes unique to control vs domatia genotypes
uniquedomcon.710 <- domcon.710[!rownames(domcon.710)%in%rownames(domcon.711),]
uniquedomcon.711 <- domcon.711[!rownames(domcon.711)%in%rownames(domcon.710),]

# Find overlap
## Shared between domatia and domatia
a <- domcon.710[rownames(domcon.710)%in%rownames(domcon.711),]

## Shared between control and domatia 710
b <- concon[rownames(concon)%in%rownames(domcon.710),]

## Shared between control and domatia 711
c <- concon[rownames(concon)%in%rownames(domcon.711),]

## Shared between domatia and domatia 710
d <- domcon.710[rownames(domcon.710)%in%rownames(domdom),]

## Shared between domatia and domatia 711
e <- domcon.711[rownames(domcon.711)%in%rownames(domdom),]

## Shared between control and domatia
f <- concon[rownames(concon)%in%rownames(domdom),]

## Shared between uniquedomcon.710 and control
g <- concon[rownames(concon)%in%rownames(uniquedomcon.710),]

## Shared between uniquedomcon.710 and domatia
h <- uniquedomcon.710[rownames(uniquedomcon.710)%in%rownames(domdom),]

## Shared between uniquedomcon.710 and control
i <- concon[rownames(concon)%in%rownames(uniquedomcon.711),]

## Shared between uniquedomcon.710 and domatia
j <- uniquedomcon.711[rownames(uniquedomcon.711)%in%rownames(domdom),]

############################# EXPLORING INTRASPECIFIC VARIATION #############################

# Get genes upregulated in X variety control compared to the other]
concon.588710.upregulated <- concon[concon$log2FoldChange<0,]
concon.588711.upregulated <- concon[concon$log2FoldChange>0,]

# Find overlap of genes upregulated in X variety domatia (compared to X variety control) with
# X variety control upregulated genes (compared to other control)
g1 <- concon.588710.upregulated[rownames(concon.588710.upregulated)%in%rownames(uniquedomcon.710),] #52/319 overlap
i1 <- concon.588711.upregulated[rownames(concon.588711.upregulated)%in%rownames(uniquedomcon.711),] #6/87 overlap

# How many of these have DECREASED expression in domatia?
uniquedomcon.710.downregulated.dom <- uniquedomcon.710[uniquedomcon.710$log2FoldChange<0,]
g1D <- concon.588710.upregulated[rownames(concon.588710.upregulated)%in%rownames(uniquedomcon.710.downregulated.dom),] #32 genes go from upregulated in leaves to downregulated in domatia

uniquedomcon.711.downregulated.dom <- uniquedomcon.711[uniquedomcon.711$log2FoldChange<0,]
i1D <- concon.588711.upregulated[rownames(concon.588711.upregulated)%in%rownames(uniquedomcon.711.downregulated.dom),] #1 gene go from upregulated in leaves to downregulated in domatia

# How many of these have INCREASED expression in domatia?
uniquedomcon.710.upregulated.dom <- uniquedomcon.710[uniquedomcon.710$log2FoldChange>0,]
g2U <- concon.588711.upregulated[rownames(concon.588711.upregulated)%in%rownames(uniquedomcon.710.upregulated.dom),] #257 genes go from downregulated in leaves to upregulated in domatia

uniquedomcon.711.upregulated.dom <- uniquedomcon.711[uniquedomcon.711$log2FoldChange>0,] 
i2U <- concon.588710.upregulated[rownames(concon.588710.upregulated)%in%rownames(uniquedomcon.711.upregulated.dom),] #75 genes go from downregulated in leaves to upregulated in domatia

# Most genes seem to be following a pattern of lowly expressed in leaf tissue (when compared to other control leaf tissue)
# to highly expressed in domatia tissue, and alternate genotype just has leaf expression=domatia expression (aka genes in g2U and i2U)

############################# GO TERM ENRICHMENT OF GENES THAT FOLLOW PATTERN IN BOTH #############################
# Load data

## Load functional annotations
#all <- read.csv("Vriparia-functional-annotations-TAIR10-mod.tsv", sep='\t')
all <- read.csv("Vriparia-functional-annotations-withparents.tsv", sep=',', na.strings=c("","NA")) # has parent already added

## Load differentially expressed genes
degenes <- unique(c(as.list(g2U$Gene), as.list(i2U$Gene)))

# Modify data

## Modify GO term file
allgo <- all[,c(2,5)]
foo<-function(x){
  inds<-allgo$gene==x
  tmp1<-unique(unlist(strsplit(allgo$all.pars[inds],"\\|"),use.names=FALSE))
  tmp1<-tmp1[nzchar(tmp1)&!is.na(tmp1)]
  data.frame("gene"=x,
             "all.pars"=paste(tmp1,collapse="|"))
}
allgo <- do.call(rbind,
                 lapply(unique(allgo$gene),foo))

allgo <- allgo %>% mutate_all(na_if,"")
allgo <- na.omit(allgo)

rn1 <- paste(allgo[,1], sep="")
gene2GO <- paste(allgo$all.pars)
names(gene2GO)<-rn1
gene2GO<-strsplit(gene2GO,"\\|")

# remove any candidate genes without GO annotation
candidate_list <- degenes
keep = candidate_list %in% allgo$gene
keep =which(keep==TRUE)
candidate_list=candidate_list[keep]

# make named factor showing which genes are of interest
bg_genes <- allgo$gene
geneList=factor(as.integer(bg_genes %in% candidate_list))
names(geneList)= bg_genes


## Set up as topGOdata with BP ontology
GODATA <- new("topGOdata",
              description = "Domatia vs Leaf in intra",
              ontology = "BP",
              allGenes = geneList,
              gene2GO = gene2GO,
              annotationFun = annFUN.gene2GO,
)

# Run GO term enrichment analysis

## Run with Fisher's exact test
resultFisher.BP <- runTest(GODATA, algorithm = "classic", statistic = "fisher")
resultFisher.BP
allGO_BP=usedGO(GODATA)
all_res_BP_intra=GenTable(GODATA, weightFisher=resultFisher.BP, orderBy='score', topNodes=length(allGO_BP))

#performing BH correction on our p values
p.adj=round(p.adjust(all_res_BP_intra$weightFisher,method="BH"),digits = 4)

# create the file with all the statistics from GO analysis
all_res_BP_intra_final=cbind(all_res_BP_intra,p.adj)
all_res_BP_intra_final=all_res_BP_intra_final[order(all_res_BP_intra_final$p.adj),]

sign_res_BP_intra <- all_res_BP_intra[all_res_BP_intra$weightFisher<0.05,]
sign_res_BP_intra_final <- all_res_BP_intra_final[all_res_BP_intra_final$p.adj<0.05,]

## Set up as topGOdata with MF ontology
GODATA <- new("topGOdata",
              description = "Domatia vs Leaf in 588intra",
              ontology = "MF",
              allGenes = geneList,
              gene2GO = gene2GO,
              annotationFun = annFUN.gene2GO,
)

# Run GO term enrichment analysis

## Run with Fisher's exact test
resultFisher.MF <- runTest(GODATA, algorithm = "classic", statistic = "fisher")
resultFisher.MF
allGO_MF=usedGO(GODATA)
all_res_MF_intra=GenTable(GODATA, weightFisher=resultFisher.MF, orderBy='score', topNodes=length(allGO_MF))
sign_res_MF_intra <- all_res_MF_intra[all_res_MF_intra$weightFisher<0.05,]

#performing BH correction on our p values
p.adj=round(p.adjust(all_res_MF_intra$weightFisher,method="BH"),digits = 4)

# create the file with all the statistics from GO analysis
all_res_MF_intra_final=cbind(all_res_MF_intra,p.adj)
all_res_MF_intra_final=all_res_MF_intra_final[order(all_res_MF_intra_final$p.adj),]

sign_res_MF_intra <- all_res_MF_intra[all_res_MF_intra$weightFisher<0.05,]
sign_res_MF_intra_final <- all_res_MF_intra_final[all_res_MF_intra_final$p.adj<0.05,]

## Set up as topGOdata with CC ontology
GODATA <- new("topGOdata",
              description = "Domatia vs Leaf in 588intra",
              ontology = "CC",
              allGenes = geneList,
              gene2GO = gene2GO,
              annotationFun = annFUN.gene2GO,
)

# Run GO term enrichment analysis

## Run with Fisher's exact test
resultFisher.CC <- runTest(GODATA, algorithm = "classic", statistic = "fisher")
resultFisher.CC
allGO_CC=usedGO(GODATA)
all_res_CC_intra=GenTable(GODATA, weightFisher=resultFisher.CC, orderBy='score', topNodes=length(allGO_CC))
sign_res_CC_intra <- all_res_CC_intra[all_res_CC_intra$weightFisher<0.05,]

#performing BH correction on our p values
p.adj=round(p.adjust(all_res_CC_intra$weightFisher,method="BH"),digits = 4)

# create the file with all the statistics from GO analysis
all_res_CC_intra_final=cbind(all_res_CC_intra,p.adj)
all_res_CC_intra_final=all_res_CC_intra_final[order(all_res_CC_intra_final$p.adj),]

sign_res_CC_intra <- all_res_CC_intra[all_res_CC_intra$weightFisher<0.05,]
sign_res_CC_intra_final <- all_res_CC_intra_final[all_res_CC_intra_final$p.adj<0.05,]

############################# GO TERM ENRICHMENT OF GENES THAT FOLLOW PATTERN IN SDG #############################
# Load data

## Load functional annotations
#all <- read.csv("Vriparia-functional-annotations-TAIR10-mod.tsv", sep='\t')
all <- read.csv("Vriparia-functional-annotations-withparents.tsv", sep=',', na.strings=c("","NA")) # has parent already added

## Load differentially expressed genes
degenes <- unique(as.list(g2U$Gene))

# Modify data

## Modify GO term file
allgo <- all[,c(2,5)]
foo<-function(x){
  inds<-allgo$gene==x
  tmp1<-unique(unlist(strsplit(allgo$all.pars[inds],"\\|"),use.names=FALSE))
  tmp1<-tmp1[nzchar(tmp1)&!is.na(tmp1)]
  data.frame("gene"=x,
             "all.pars"=paste(tmp1,collapse="|"))
}
allgo <- do.call(rbind,
                 lapply(unique(allgo$gene),foo))

allgo <- allgo %>% mutate_all(na_if,"")
allgo <- na.omit(allgo)

rn1 <- paste(allgo[,1], sep="")
gene2GO <- paste(allgo$all.pars)
names(gene2GO)<-rn1
gene2GO<-strsplit(gene2GO,"\\|")

# remove any candidate genes without GO annotation
candidate_list <- degenes
keep = candidate_list %in% allgo$gene
keep =which(keep==TRUE)
candidate_list=candidate_list[keep]

# make named factor showing which genes are of interest
bg_genes <- allgo$gene
geneList=factor(as.integer(bg_genes %in% candidate_list))
names(geneList)= bg_genes


## Set up as topGOdata with BP ontology
GODATA <- new("topGOdata",
              description = "Domatia vs Leaf in intra.710",
              ontology = "BP",
              allGenes = geneList,
              gene2GO = gene2GO,
              annotationFun = annFUN.gene2GO,
)

# Run GO term enrichment analysis

## Run with Fisher's exact test
resultFisher.BP <- runTest(GODATA, algorithm = "classic", statistic = "fisher")
resultFisher.BP
allGO_BP=usedGO(GODATA)
all_res_BP_intra.710=GenTable(GODATA, weightFisher=resultFisher.BP, orderBy='score', topNodes=length(allGO_BP))

#performing BH correction on our p values
p.adj=round(p.adjust(all_res_BP_intra.710$weightFisher,method="BH"),digits = 4)

# create the file with all the statistics from GO analysis
all_res_BP_intra.710_final=cbind(all_res_BP_intra.710,p.adj)
all_res_BP_intra.710_final=all_res_BP_intra.710_final[order(all_res_BP_intra.710_final$p.adj),]

sign_res_BP_intra.710 <- all_res_BP_intra.710[all_res_BP_intra.710$weightFisher<0.05,]
sign_res_BP_intra.710_final <- all_res_BP_intra.710_final[all_res_BP_intra.710_final$p.adj<0.05,]

## Set up as topGOdata with MF ontology
GODATA <- new("topGOdata",
              description = "Domatia vs Leaf in 588intra.710",
              ontology = "MF",
              allGenes = geneList,
              gene2GO = gene2GO,
              annotationFun = annFUN.gene2GO,
)

# Run GO term enrichment analysis

## Run with Fisher's exact test
resultFisher.MF <- runTest(GODATA, algorithm = "classic", statistic = "fisher")
resultFisher.MF
allGO_MF=usedGO(GODATA)
all_res_MF_intra.710=GenTable(GODATA, weightFisher=resultFisher.MF, orderBy='score', topNodes=length(allGO_MF))
sign_res_MF_intra.710 <- all_res_MF_intra.710[all_res_MF_intra.710$weightFisher<0.05,]

#performing BH correction on our p values
p.adj=round(p.adjust(all_res_MF_intra.710$weightFisher,method="BH"),digits = 4)

# create the file with all the statistics from GO analysis
all_res_MF_intra.710_final=cbind(all_res_MF_intra.710,p.adj)
all_res_MF_intra.710_final=all_res_MF_intra.710_final[order(all_res_MF_intra.710_final$p.adj),]

sign_res_MF_intra.710 <- all_res_MF_intra.710[all_res_MF_intra.710$weightFisher<0.05,]
sign_res_MF_intra.710_final <- all_res_MF_intra.710_final[all_res_MF_intra.710_final$p.adj<0.05,]

## Set up as topGOdata with CC ontology
GODATA <- new("topGOdata",
              description = "Domatia vs Leaf in 588intra.710",
              ontology = "CC",
              allGenes = geneList,
              gene2GO = gene2GO,
              annotationFun = annFUN.gene2GO,
)

# Run GO term enrichment analysis

## Run with Fisher's exact test
resultFisher.CC <- runTest(GODATA, algorithm = "classic", statistic = "fisher")
resultFisher.CC
allGO_CC=usedGO(GODATA)
all_res_CC_intra.710=GenTable(GODATA, weightFisher=resultFisher.CC, orderBy='score', topNodes=length(allGO_CC))
sign_res_CC_intra.710 <- all_res_CC_intra.710[all_res_CC_intra.710$weightFisher<0.05,]

#performing BH correction on our p values
p.adj=round(p.adjust(all_res_CC_intra.710$weightFisher,method="BH"),digits = 4)

# create the file with all the statistics from GO analysis
all_res_CC_intra.710_final=cbind(all_res_CC_intra.710,p.adj)
all_res_CC_intra.710_final=all_res_CC_intra.710_final[order(all_res_CC_intra.710_final$p.adj),]

sign_res_CC_intra.710 <- all_res_CC_intra.710[all_res_CC_intra.710$weightFisher<0.05,]
sign_res_CC_intra.710_final <- all_res_CC_intra.710_final[all_res_CC_intra.710_final$p.adj<0.05,]

# Plot BP results
# Look at overlap for all significantly enriched GO terms
sign_res_BP_intra.710_final$Genotype <- c("SDG")
sign_res_BP_intra.710_final$Log.Fold.Enrichment <- log(sign_res_BP_intra.710_final$Significant / sign_res_BP_intra.710_final$Expected)
sign_res_BP_intra.710_final$fullbeans <-  paste0(sign_res_BP_intra.710_final$Term, sep = " (", sign_res_BP_intra.710_final$GO.ID, sep = ")")

ggplot(sign_res_BP_intra.710_final, (aes(x=Genotype, y=fullbeans, color = as.numeric(Log.Fold.Enrichment), size=as.numeric(Significant)))) + 
  geom_point() +
  scale_color_gradient(low = "orange", high = "blue", name = "Log Fold Enrichment") +
  scale_y_discrete(limits=rev) +
  scale_x_discrete(limits=rev) +
  theme_classic() +
  ylab("GO Term (Biological Process)") +
  scale_size_continuous(name = "DEGs") +
  theme(axis.text.x = element_text(size=11, face="bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=11, face="bold"))

############################# GO TERM ENRICHMENT OF GENES THAT FOLLOW PATTERN IN LDG #############################
# Load data

## Load functional annotations
#all <- read.csv("Vriparia-functional-annotations-TAIR10-mod.tsv", sep='\t')
all <- read.csv("Vriparia-functional-annotations-withparents.tsv", sep=',', na.strings=c("","NA")) # has parent already added

## Load differentially expressed genes
degenes <- unique(as.list(i2U$Gene))

# Modify data

## Modify GO term file
allgo <- all[,c(2,5)]
foo<-function(x){
  inds<-allgo$gene==x
  tmp1<-unique(unlist(strsplit(allgo$all.pars[inds],"\\|"),use.names=FALSE))
  tmp1<-tmp1[nzchar(tmp1)&!is.na(tmp1)]
  data.frame("gene"=x,
             "all.pars"=paste(tmp1,collapse="|"))
}
allgo <- do.call(rbind,
                 lapply(unique(allgo$gene),foo))

allgo <- allgo %>% mutate_all(na_if,"")
allgo <- na.omit(allgo)

rn1 <- paste(allgo[,1], sep="")
gene2GO <- paste(allgo$all.pars)
names(gene2GO)<-rn1
gene2GO<-strsplit(gene2GO,"\\|")

# remove any candidate genes without GO annotation
candidate_list <- degenes
keep = candidate_list %in% allgo$gene
keep =which(keep==TRUE)
candidate_list=candidate_list[keep]

# make named factor showing which genes are of interest
bg_genes <- allgo$gene
geneList=factor(as.integer(bg_genes %in% candidate_list))
names(geneList)= bg_genes


## Set up as topGOdata with BP ontology
GODATA <- new("topGOdata",
              description = "Domatia vs Leaf in intra.711",
              ontology = "BP",
              allGenes = geneList,
              gene2GO = gene2GO,
              annotationFun = annFUN.gene2GO,
)

# Run GO term enrichment analysis

## Run with Fisher's exact test
resultFisher.BP <- runTest(GODATA, algorithm = "classic", statistic = "fisher")
resultFisher.BP
allGO_BP=usedGO(GODATA)
all_res_BP_intra.711=GenTable(GODATA, weightFisher=resultFisher.BP, orderBy='score', topNodes=length(allGO_BP))

#performing BH correction on our p values
p.adj=round(p.adjust(all_res_BP_intra.711$weightFisher,method="BH"),digits = 4)

# create the file with all the statistics from GO analysis
all_res_BP_intra.711_final=cbind(all_res_BP_intra.711,p.adj)
all_res_BP_intra.711_final=all_res_BP_intra.711_final[order(all_res_BP_intra.711_final$p.adj),]

sign_res_BP_intra.711 <- all_res_BP_intra.711[all_res_BP_intra.711$weightFisher<0.05,]
sign_res_BP_intra.711_final <- all_res_BP_intra.711_final[all_res_BP_intra.711_final$p.adj<0.05,]

## Set up as topGOdata with MF ontology
GODATA <- new("topGOdata",
              description = "Domatia vs Leaf in 588intra.711",
              ontology = "MF",
              allGenes = geneList,
              gene2GO = gene2GO,
              annotationFun = annFUN.gene2GO,
)

# Run GO term enrichment analysis

## Run with Fisher's exact test
resultFisher.MF <- runTest(GODATA, algorithm = "classic", statistic = "fisher")
resultFisher.MF
allGO_MF=usedGO(GODATA)
all_res_MF_intra.711=GenTable(GODATA, weightFisher=resultFisher.MF, orderBy='score', topNodes=length(allGO_MF))
sign_res_MF_intra.711 <- all_res_MF_intra.711[all_res_MF_intra.711$weightFisher<0.05,]

#performing BH correction on our p values
p.adj=round(p.adjust(all_res_MF_intra.711$weightFisher,method="BH"),digits = 4)

# create the file with all the statistics from GO analysis
all_res_MF_intra.711_final=cbind(all_res_MF_intra.711,p.adj)
all_res_MF_intra.711_final=all_res_MF_intra.711_final[order(all_res_MF_intra.711_final$p.adj),]

sign_res_MF_intra.711 <- all_res_MF_intra.711[all_res_MF_intra.711$weightFisher<0.05,]
sign_res_MF_intra.711_final <- all_res_MF_intra.711_final[all_res_MF_intra.711_final$p.adj<0.05,]

## Set up as topGOdata with CC ontology
GODATA <- new("topGOdata",
              description = "Domatia vs Leaf in 588intra.711",
              ontology = "CC",
              allGenes = geneList,
              gene2GO = gene2GO,
              annotationFun = annFUN.gene2GO,
)

# Run GO term enrichment analysis

## Run with Fisher's exact test
resultFisher.CC <- runTest(GODATA, algorithm = "classic", statistic = "fisher")
resultFisher.CC
allGO_CC=usedGO(GODATA)
all_res_CC_intra.711=GenTable(GODATA, weightFisher=resultFisher.CC, orderBy='score', topNodes=length(allGO_CC))
sign_res_CC_intra.711 <- all_res_CC_intra.711[all_res_CC_intra.711$weightFisher<0.05,]

#performing BH correction on our p values
p.adj=round(p.adjust(all_res_CC_intra.711$weightFisher,method="BH"),digits = 4)

# create the file with all the statistics from GO analysis
all_res_CC_intra.711_final=cbind(all_res_CC_intra.711,p.adj)
all_res_CC_intra.711_final=all_res_CC_intra.711_final[order(all_res_CC_intra.711_final$p.adj),]

sign_res_CC_intra.711 <- all_res_CC_intra.711[all_res_CC_intra.711$weightFisher<0.05,]
sign_res_CC_intra.711_final <- all_res_CC_intra.711_final[all_res_CC_intra.711_final$p.adj<0.05,]