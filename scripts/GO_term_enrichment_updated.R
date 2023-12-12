# Set working directory and load necessary packages
setwd("C:/Users/rittere5/OneDrive - Michigan State University/Vitis-domatia/")
setwd("C:/Users/elean/OneDrive - Michigan State University/Vitis-domatia/")

library("topGO")
library(dplyr)
library(tidyr)
library(ggplot2)

# NOTE: a version of the old code that contains the code for NOT using parents (as well as the 
# code that was used to get the parents) back in is at the bottom, saved for a rainy day

######################## GO TERM ENRICHMENT WITH DOMATIA V LEAF IN 588710 ######################## 

# Load data

## Load functional annotations
#all <- read.csv("Vriparia-functional-annotations-TAIR10-mod.tsv", sep='\t')
all <- read.csv("Vriparia-functional-annotations-withparents.tsv", sep=',', na.strings=c("","NA")) # has parent already added

## Load differentially expressed genes
degenes <- read.csv("differentially_expressed_genes/DE_genes_Domatia_V_Leaf_588710.csv")
degenes <- degenes[degenes$log2FoldChange>0 ,] #Keep only upregulated genes

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
candidate_list <- degenes$Gene
keep = candidate_list %in% allgo$gene
keep =which(keep==TRUE)
candidate_list=candidate_list[keep]

# make named factor showing which genes are of interest
bg_genes <- allgo$gene
geneList=factor(as.integer(bg_genes %in% candidate_list))
names(geneList)= bg_genes


## Set up as topGOdata with BP ontology
GODATA <- new("topGOdata",
              description = "Domatia vs Leaf in 588710",
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
all_res_BP_710=GenTable(GODATA, weightFisher=resultFisher.BP, orderBy='score', topNodes=length(allGO_BP))

#performing BH correction on our p values
p.adj=round(p.adjust(all_res_BP_710$weightFisher,method="BH"),digits = 4)

# create the file with all the statistics from GO analysis
all_res_BP_710_final=cbind(all_res_BP_710,p.adj)
all_res_BP_710_final=all_res_BP_710_final[order(all_res_BP_710_final$p.adj),]

sign_res_BP_710 <- all_res_BP_710[all_res_BP_710$weightFisher<0.05,]
sign_res_BP_710_final <- all_res_BP_710_final[all_res_BP_710_final$p.adj<0.05,]

#write.csv(all_res_BP_710_final, "GO-term-enrichment/Upregulated-genes-588710-Con-vs-Dom-GO-terms-BP.csv", row.names = FALSE)

## Set up as topGOdata with MF ontology
GODATA <- new("topGOdata",
              description = "Domatia vs Leaf in 588710",
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
all_res_MF_710=GenTable(GODATA, weightFisher=resultFisher.MF, orderBy='score', topNodes=length(allGO_MF))
sign_res_MF_710 <- all_res_MF_710[all_res_MF_710$weightFisher<0.05,]

#performing BH correction on our p values
p.adj=round(p.adjust(all_res_MF_710$weightFisher,method="BH"),digits = 4)

# create the file with all the statistics from GO analysis
all_res_MF_710_final=cbind(all_res_MF_710,p.adj)
all_res_MF_710_final=all_res_MF_710_final[order(all_res_MF_710_final$p.adj),]

sign_res_MF_710 <- all_res_MF_710[all_res_MF_710$weightFisher<0.05,]
sign_res_MF_710_final <- all_res_MF_710_final[all_res_MF_710_final$p.adj<0.05,]

#write.csv(all_res_MF_710_final, "GO-term-enrichment/Upregulated-genes-588710-Con-vs-Dom-GO-terms-MF.csv", row.names = FALSE)

## Set up as topGOdata with CC ontology
GODATA <- new("topGOdata",
              description = "Domatia vs Leaf in 588710",
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
all_res_CC_710=GenTable(GODATA, weightFisher=resultFisher.CC, orderBy='score', topNodes=length(allGO_CC))
sign_res_CC_710 <- all_res_CC_710[all_res_CC_710$weightFisher<0.05,]

#performing BH correction on our p values
p.adj=round(p.adjust(all_res_CC_710$weightFisher,method="BH"),digits = 4)

# create the file with all the statistics from GO analysis
all_res_CC_710_final=cbind(all_res_CC_710,p.adj)
all_res_CC_710_final=all_res_CC_710_final[order(all_res_CC_710_final$p.adj),]

sign_res_CC_710 <- all_res_CC_710[all_res_CC_710$weightFisher<0.05,]
sign_res_CC_710_final <- all_res_CC_710_final[all_res_CC_710_final$p.adj<0.05,]

#write.csv(all_res_CC_710_final, "GO-term-enrichment/Upregulated-genes-588710-Con-vs-Dom-GO-terms-CC.csv", row.names = FALSE)

######################## GO TERM ENRICHMENT WITH DOMATIA V LEAF IN 588711 ######################## 

# Load data

## Load functional annotations
#all <- read.csv("Vriparia-functional-annotations-TAIR10-mod.tsv", sep='\t')
all <- read.csv("Vriparia-functional-annotations-withparents.tsv", sep=',', na.strings=c("","NA")) # has parent already added

## Load differentially expressed genes
degenes <- read.csv("differentially_expressed_genes/DE_genes_Domatia_V_Leaf_588711.csv")
degenes <- degenes[degenes$log2FoldChange>0 ,] #Keep only upregulated genes

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
candidate_list <- degenes$Gene
keep = candidate_list %in% allgo$gene
keep =which(keep==TRUE)
candidate_list=candidate_list[keep]

# make named factor showing which genes are of interest
bg_genes <- allgo$gene
geneList=factor(as.integer(bg_genes %in% candidate_list))
names(geneList)= bg_genes


## Set up as topGOdata with BP ontology
GODATA <- new("topGOdata",
              description = "Domatia vs Leaf in 588711",
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
all_res_BP_711=GenTable(GODATA, weightFisher=resultFisher.BP, orderBy='score', topNodes=length(allGO_BP))
sign_res_BP_711 <- all_res_BP_711[all_res_BP_711$weightFisher<0.05,]

#performing BH correction on our p values
p.adj=round(p.adjust(all_res_BP_711$weightFisher,method="BH"),digits = 4)

# create the file with all the statistics from GO analysis
all_res_BP_711_final=cbind(all_res_BP_711,p.adj)
all_res_BP_711_final=all_res_BP_711_final[order(all_res_BP_711_final$p.adj),]

sign_res_BP_711 <- all_res_BP_711[all_res_BP_711$weightFisher<0.05,]
sign_res_BP_711_final <- all_res_BP_711_final[all_res_BP_711_final$p.adj<0.05,]

#write.csv(all_res_BP_711_final, "GO-term-enrichment/Upregulated-genes-588711-Con-vs-Dom-GO-terms-BP.csv", row.names = FALSE)

## Set up as topGOdata with MF ontology
GODATA <- new("topGOdata",
              description = "Domatia vs Leaf in 588711",
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
all_res_MF_711=GenTable(GODATA, weightFisher=resultFisher.MF, orderBy='score', topNodes=length(allGO_MF))
sign_res_MF_711 <- all_res_MF_711[all_res_MF_711$weightFisher<0.05,]

#performing BH correction on our p values
p.adj=round(p.adjust(all_res_MF_711$weightFisher,method="BH"),digits = 4)

# create the file with all the statistics from GO analysis
all_res_MF_711_final=cbind(all_res_MF_711,p.adj)
all_res_MF_711_final=all_res_MF_711_final[order(all_res_MF_711_final$p.adj),]

sign_res_MF_711 <- all_res_MF_711[all_res_MF_711$weightFisher<0.05,]
sign_res_MF_711_final <- all_res_MF_711_final[all_res_MF_711_final$p.adj<0.05,]

#write.csv(all_res_MF_711_final, "GO-term-enrichment/Upregulated-genes-588711-Con-vs-Dom-GO-terms-MF.csv", row.names = FALSE)

## Set up as topGOdata with CC ontology
GODATA <- new("topGOdata",
              description = "Domatia vs Leaf in 588711",
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
all_res_CC_711=GenTable(GODATA, weightFisher=resultFisher.CC, orderBy='score', topNodes=length(allGO_CC))
sign_res_CC_711 <- all_res_CC_711[all_res_CC_711$weightFisher<0.05,]


#performing BH correction on our p values
p.adj=round(p.adjust(all_res_CC_711$weightFisher,method="BH"),digits = 4)

# create the file with all the statistics from GO analysis
all_res_CC_711_final=cbind(all_res_CC_711,p.adj)
all_res_CC_711_final=all_res_CC_711_final[order(all_res_CC_711_final$p.adj),]

sign_res_CC_711 <- all_res_CC_711[all_res_CC_711$weightFisher<0.05,]
sign_res_CC_711_final <- all_res_CC_711_final[all_res_CC_711_final$p.adj<0.05,]

#write.csv(all_res_CC_711_final, "GO-term-enrichment/Upregulated-genes-588711-Con-vs-Dom-GO-terms-CC.csv", row.names = FALSE)

######################## OVERLAP BETWEEN ENRICHED GO TERMS ########################
# BP
bp.overlap <- merge(sign_res_BP_710_final, sign_res_BP_711_final, by="GO.ID")
writeClipboard(paste0(as.vector(bp.overlap$GO.ID)))

# MF
mf.overlap <- merge(sign_res_MF_710_final, sign_res_MF_711_final, by="GO.ID")
writeClipboard(paste0(as.vector(mf.overlap$GO.ID)))

# CC
cc.overlap <- merge(sign_res_CC_710_final, sign_res_CC_711_final, by="GO.ID")
writeClipboard(paste0(as.vector(cc.overlap$GO.ID)))

######################## MAKE PLOT OF OVERLAPPING BP TERMS ########################
all_res_BP_710_final <- read.csv("GO-term-enrichment/Upregulated-genes-588710-Con-vs-Dom-GO-terms-BP.csv")
all_res_BP_711_final <- read.csv("GO-term-enrichment/Upregulated-genes-588711-Con-vs-Dom-GO-terms-BP.csv")

sign_res_BP_710_final <- all_res_BP_710_final[all_res_BP_710_final$p.adj<0.05 ,]
sign_res_BP_711_final <- all_res_BP_711_final[all_res_BP_711_final$p.adj<0.05 ,]

# Look at overlap for all significantly enriched GO terms
sign_res_BP_710_final$Genotype <- c("SDG")
sign_res_BP_711_final$Genotype <- c("LDG")
tempa <- sign_res_BP_710_final[sign_res_BP_710_final$GO.ID %in% sign_res_BP_711_final$GO.ID ,]
tempb <- sign_res_BP_711_final[sign_res_BP_711_final$GO.ID %in% sign_res_BP_710_final$GO.ID ,]
sign.overlap <- unique(rbind(tempa, tempb))
sign.overlap$Log.Fold.Enrichment <- log(sign.overlap$Significant / sign.overlap$Expected)
sign.overlap <- data.frame(lapply(sign.overlap,
                                  function(x) gsub("anatomical structure formation involved ...",
                                                   "anatomical structure formation involved in morphogenesis", x)))
sign.overlap$fullbeans <-  paste0(sign.overlap$Term, sep = " (", sign.overlap$GO.ID, sep = ")")

ggplot(sign.overlap, (aes(x=Genotype, y=fullbeans, color = as.numeric(Log.Fold.Enrichment), size=as.numeric(Significant)))) + 
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

######################## MAKE PLOT OF OVERLAPPING MF TERMS ########################
all_res_MF_710_final <- read.csv("GO-term-enrichment/Upregulated-genes-588710-Con-vs-Dom-GO-terms-MF.csv")
all_res_MF_711_final <- read.csv("GO-term-enrichment/Upregulated-genes-588711-Con-vs-Dom-GO-terms-MF.csv")

sign_res_MF_710_final <- all_res_MF_710_final[all_res_MF_710_final$p.adj<0.05 ,]
sign_res_MF_711_final <- all_res_MF_711_final[all_res_MF_711_final$p.adj<0.05 ,]

# Look at overlap for all significantly enriched GO terms
sign_res_MF_710_final$Genotype <- c("SDG")
sign_res_MF_711_final$Genotype <- c("LDG")
tempa <- sign_res_MF_710_final[sign_res_MF_710_final$GO.ID %in% sign_res_MF_711_final$GO.ID ,]
tempb <- sign_res_MF_711_final[sign_res_MF_711_final$GO.ID %in% sign_res_MF_710_final$GO.ID ,]
sign.overlap <- unique(rbind(tempa, tempb))
sign.overlap$Log.Fold.Enrichment <- log(sign.overlap$Significant / sign.overlap$Expected)
sign.overlap <- data.frame(lapply(sign.overlap,
                                  function(x) gsub("amino acid transmembrane transporter act...",
                                                   "amino acid transmembrane transporter activity", x)))
sign.overlap <- data.frame(lapply(sign.overlap,
                                  function(x) gsub("carboxylic acid transmembrane transporte...",
                                                   "carboxylic acid transmembrane transporter activity", x)))
sign.overlap <- data.frame(lapply(sign.overlap,
                                  function(x) gsub("inorganic molecular entity transmembrane...",
                                                   "inorganic molecular entity transmembrane transporter activity", x)))
sign.overlap <- data.frame(lapply(sign.overlap,
                                  function(x) gsub("organic acid transmembrane transporter a...",
                                                   "organic acid transmembrane transporter activity", x)))
sign.overlap <- data.frame(lapply(sign.overlap,
                                  function(x) gsub("sequence-specific double-stranded DNA bi...",
                                                   "sequence-specific double-stranded DNA binding", x)))
sign.overlap <- data.frame(lapply(sign.overlap,
                                  function(x) gsub("transcription regulatory region sequence...",
                                                   "transcription regulatory region sequence-specific DNA binding", x)))
sign.overlap$fullbeans <-  paste0(sign.overlap$Term, sep = " (", sign.overlap$GO.ID, sep = ")")

ggplot(sign.overlap, (aes(x=Genotype, y=fullbeans, color = as.numeric(Log.Fold.Enrichment), size=as.numeric(Significant)))) + 
  geom_point() +
  scale_color_gradient(low = "orange", high = "blue", name = "Log Fold Enrichment") +
  scale_y_discrete(limits=rev) +
  scale_x_discrete(limits=rev) +
  theme_classic() +
  ylab("GO Term (Molecular Function)") +
  scale_size_continuous(name = "DEGs") +
  theme(axis.text.x = element_text(size=11, face="bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=11, face="bold"))

######################## MAKE PLOT OF OVERLAPPING CC TERMS ########################
all_res_CC_710_final <- read.csv("GO-term-enrichment/Upregulated-genes-588710-Con-vs-Dom-GO-terms-CC.csv")
all_res_CC_711_final <- read.csv("GO-term-enrichment/Upregulated-genes-588711-Con-vs-Dom-GO-terms-CC.csv")

sign_res_CC_710_final <- all_res_CC_710_final[all_res_CC_710_final$p.adj<0.05 ,]
sign_res_CC_711_final <- all_res_CC_711_final[all_res_CC_711_final$p.adj<0.05 ,]

# Look at overlap for all significantly enriched GO terms
sign_res_CC_710_final$Genotype <- c("SDG")
sign_res_CC_711_final$Genotype <- c("LDG")
tempa <- sign_res_CC_710_final[sign_res_CC_710_final$GO.ID %in% sign_res_CC_711_final$GO.ID ,]
tempb <- sign_res_CC_711_final[sign_res_CC_711_final$GO.ID %in% sign_res_CC_710_final$GO.ID ,]
sign.overlap <- unique(rbind(tempa, tempb))
sign.overlap$Log.Fold.Enrichment <- log(sign.overlap$Significant / sign.overlap$Expected)
sign.overlap$fullbeans <-  paste0(sign.overlap$Term, sep = " (", sign.overlap$GO.ID, sep = ")")

ggplot(sign.overlap, (aes(x=Genotype, y=fullbeans, color = as.numeric(Log.Fold.Enrichment), size=as.numeric(Significant)))) + 
  geom_point() +
  scale_color_gradient(low = "orange", high = "blue", name = "Log Fold Enrichment") +
  scale_y_discrete(limits=rev) +
  scale_x_discrete(limits=rev) +
  theme_classic() +
  ylab("GO Term (Cellular Component)") +
  scale_size_continuous(name = "DEGs") +
  theme(axis.text.x = element_text(size=11, face="bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=11, face="bold"))

######################## GO TERM ENRICHMENT WITH DOMATIA V DOMATIA ######################## 

# Load data

## Load functional annotations
#all <- read.csv("Vriparia-functional-annotations-TAIR10-mod.tsv", sep='\t')
all <- read.csv("Vriparia-functional-annotations-withparents.tsv", sep=',', na.strings=c("","NA")) # has parent already added

## Load differentially expressed genes
degenes <- read.csv("differentially_expressed_genes/DE_genes_Domatia_710_v_711.csv")
degenes <- degenes[degenes$log2FoldChange>0 ,] #Keep only upregulated genes

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
candidate_list <- degenes$Gene
keep = candidate_list %in% allgo$gene
keep =which(keep==TRUE)
candidate_list=candidate_list[keep]

# make named factor showing which genes are of interest
bg_genes <- allgo$gene
geneList=factor(as.integer(bg_genes %in% candidate_list))
names(geneList)= bg_genes


## Set up as topGOdata with BP ontology
GODATA <- new("topGOdata",
              description = "Domatia vs Leaf in 588dom",
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
all_res_BP_dom=GenTable(GODATA, weightFisher=resultFisher.BP, orderBy='score', topNodes=length(allGO_BP))
sign_res_BP_dom <- all_res_BP_dom[all_res_BP_dom$weightFisher<0.05,]

#performing BH correction on our p values
p.adj=round(p.adjust(all_res_BP_dom$weightFisher,method="BH"),digits = 4)

# create the file with all the statistics from GO analysis
all_res_BP_dom_final=cbind(all_res_BP_dom,p.adj)
all_res_BP_dom_final=all_res_BP_dom_final[order(all_res_BP_dom_final$p.adj),]

sign_res_BP_dom <- all_res_BP_dom[all_res_BP_dom$weightFisher<0.05,]
sign_res_BP_dom_final <- all_res_BP_dom_final[all_res_BP_dom_final$p.adj<0.05,]

## Set up as topGOdata with MF ontology
GODATA <- new("topGOdata",
              description = "Domatia vs Leaf in 588dom",
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
all_res_MF_dom=GenTable(GODATA, weightFisher=resultFisher.MF, orderBy='score', topNodes=length(allGO_MF))
sign_res_MF_dom <- all_res_MF_dom[all_res_MF_dom$weightFisher<0.05,]

#performing BH correction on our p values
p.adj=round(p.adjust(all_res_MF_dom$weightFisher,method="BH"),digits = 4)

# create the file with all the statistics from GO analysis
all_res_MF_dom_final=cbind(all_res_MF_dom,p.adj)
all_res_MF_dom_final=all_res_MF_dom_final[order(all_res_MF_dom_final$p.adj),]

sign_res_MF_dom <- all_res_MF_dom[all_res_MF_dom$weightFisher<0.05,]
sign_res_MF_dom_final <- all_res_MF_dom_final[all_res_MF_dom_final$p.adj<0.05,]

## Set up as topGOdata with CC ontology
GODATA <- new("topGOdata",
              description = "Domatia vs Leaf in 588dom",
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
all_res_CC_dom=GenTable(GODATA, weightFisher=resultFisher.CC, orderBy='score', topNodes=length(allGO_CC))
sign_res_CC_dom <- all_res_CC_dom[all_res_CC_dom$weightFisher<0.05,]


#performing BH correction on our p values
p.adj=round(p.adjust(all_res_CC_dom$weightFisher,method="BH"),digits = 4)

# create the file with all the statistics from GO analysis
all_res_CC_dom_final=cbind(all_res_CC_dom,p.adj)
all_res_CC_dom_final=all_res_CC_dom_final[order(all_res_CC_dom_final$p.adj),]

sign_res_CC_dom <- all_res_CC_dom[all_res_CC_dom$weightFisher<0.05,]
sign_res_CC_dom_final <- all_res_CC_dom_final[all_res_CC_dom_final$p.adj<0.05,]

######################## GO TERM ENRICHMENT WITH CONTROL V CONTROL ######################## 

# Load data

## Load functional annotations
#all <- read.csv("Vriparia-functional-annotations-TAIR10-mod.tsv", sep='\t')
all <- read.csv("Vriparia-functional-annotations-withparents.tsv", sep=',', na.strings=c("","NA")) # has parent already added

## Load differentially expressed genes
degenes <- read.csv("differentially_expressed_genes/DE_genes_Control_710_v_711.csv")
degenes <- degenes[degenes$log2FoldChange>0 ,] #Keep only upregulated genes

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
#gene2GO <- paste(allgo$Arabidopsis_GO_terms, allgo$PFAM_GO_terms, sep="|")
#gene2GO = allgo[,-1]
gene2GO <- paste(allgo$all.pars)
names(gene2GO)<-rn1
gene2GO<-strsplit(gene2GO,"\\|")

# remove any candidate genes without GO annotation
candidate_list <- degenes$Gene
keep = candidate_list %in% allgo$gene
keep =which(keep==TRUE)
candidate_list=candidate_list[keep]

# make named factor showing which genes are of interest
bg_genes <- allgo$gene
geneList=factor(as.integer(bg_genes %in% candidate_list))
names(geneList)= bg_genes


## Set up as topGOdata with BP ontology
GODATA <- new("topGOdata",
              description = "conatia vs Leaf in 588con",
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
all_res_BP_con=GenTable(GODATA, weightFisher=resultFisher.BP, orderBy='score', topNodes=length(allGO_BP))
sign_res_BP_con <- all_res_BP_con[all_res_BP_con$weightFisher<0.05,]

#performing BH correction on our p values
p.adj=round(p.adjust(all_res_BP_con$weightFisher,method="BH"),digits = 4)

# create the file with all the statistics from GO analysis
all_res_BP_con_final=cbind(all_res_BP_con,p.adj)
all_res_BP_con_final=all_res_BP_con_final[order(all_res_BP_con_final$p.adj),]

sign_res_BP_con <- all_res_BP_con[all_res_BP_con$weightFisher<0.05,]
sign_res_BP_con_final <- all_res_BP_con_final[all_res_BP_con_final$p.adj<0.05,]

## Set up as topGOdata with MF ontology
GODATA <- new("topGOdata",
              description = "conatia vs Leaf in 588con",
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
all_res_MF_con=GenTable(GODATA, weightFisher=resultFisher.MF, orderBy='score', topNodes=length(allGO_MF))
sign_res_MF_con <- all_res_MF_con[all_res_MF_con$weightFisher<0.05,]

#performing BH correction on our p values
p.adj=round(p.adjust(all_res_MF_con$weightFisher,method="BH"),digits = 4)

# create the file with all the statistics from GO analysis
all_res_MF_con_final=cbind(all_res_MF_con,p.adj)
all_res_MF_con_final=all_res_MF_con_final[order(all_res_MF_con_final$p.adj),]

sign_res_MF_con <- all_res_MF_con[all_res_MF_con$weightFisher<0.05,]
sign_res_MF_con_final <- all_res_MF_con_final[all_res_MF_con_final$p.adj<0.05,]

## Set up as topGOdata with CC ontology
GODATA <- new("topGOdata",
              description = "conatia vs Leaf in 588con",
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
all_res_CC_con=GenTable(GODATA, weightFisher=resultFisher.CC, orderBy='score', topNodes=length(allGO_CC))
sign_res_CC_con <- all_res_CC_con[all_res_CC_con$weightFisher<0.05,]


#performing BH correction on our p values
p.adj=round(p.adjust(all_res_CC_con$weightFisher,method="BH"),digits = 4)

# create the file with all the statistics from GO analysis
all_res_CC_con_final=cbind(all_res_CC_con,p.adj)
all_res_CC_con_final=all_res_CC_con_final[order(all_res_CC_con_final$p.adj),]

sign_res_CC_con <- all_res_CC_con[all_res_CC_con$weightFisher<0.05,]
sign_res_CC_con_final <- all_res_CC_con_final[all_res_CC_con_final$p.adj<0.05,]

######################## UNUSED CODE UNLESS THINGS CHANGE ######################## 

# Load data

## Load functional annotations
#all <- read.csv("Vriparia-functional-annotations-TAIR10-mod.tsv", sep='\t')
all <- read.csv("Vriparia-functional-annotations-withparents.tsv", sep=',', na.strings=c("","NA")) # has parent already added

## Add gene data to all if needed
# dict <- read.csv("cds-to-gene-V2.tsv", sep = '\t')
# colnames(dict) <- c("Transcript", "gene")
# dict$Transcript <- gsub("\\.[0-9]", "", dict$Transcript)
# all <- merge(all, dict, by="Transcript", all = TRUE)
# all <- all[,c(1, 8, 3:7)]
# all <- all %>% drop_na(gene)
# all <- all[!duplicated(all$gene),] # Get rid of duplicates - all rows are identical (transcripts are sometimes different, that is all)
# 
# ## Proving that duplicated rows are identical for all gene names and GO terms as a sanity check
# 
# test <- all
# check<-do.call(rbind,tapply(1:nrow(test),
#                             test$gene,
#                             function(ii) lengths(lapply(test,function(jj) unique(jj[ii])))))
# apply(check,2,function(ii) any(ii>2))

## Load differentially expressed genes
degenes <- read.csv("differentially_expressed_genes/DE_genes_Control_710_v_711.csv")
degenes <- degenes[degenes$log2FoldChange>0 ,] #Keep only upregulated genes

# Modify data

## Modify GO term file
#allgo <- all[,c(2,4,6)]
allgo <- all[,c(2,5)]
foo<-function(x){
  inds<-allgo$gene==x
  tmp1<-unique(unlist(strsplit(allgo$all.pars[inds],"\\|"),use.names=FALSE))
  #tmp1<-unique(unlist(strsplit(allgo$Arabidopsis_GO_terms[inds],"\\|"),use.names=FALSE))
  tmp1<-tmp1[nzchar(tmp1)&!is.na(tmp1)]
  #tmp2<-unique(unlist(strsplit(allgo$PFAM_GO_terms[inds],"\\|"),use.names=FALSE))
  #tmp2<-tmp2[nzchar(tmp2)&!is.na(tmp2)]
  data.frame("gene"=x,
             "all.pars"=paste(tmp1,collapse="|"))
  #"Arabidopsis_GO_terms"=paste(tmp1,collapse="|"),
  #"PFAM_GO_terms"=paste(tmp2,collapse="|"))
}
allgo <- do.call(rbind,
                 lapply(unique(allgo$gene),foo))

allgo <- allgo %>% mutate_all(na_if,"")
allgo <- na.omit(allgo)

## Add back in parent GO terms

# #parse obo file in lookup list
# tmp<-readLines("go-basic.obo")
# tmp<-tmp[nchar(tmp)>0]
# tmp<-tmp[do.call(`:`,as.list(match(c("[Term]","[Typedef]"),tmp)-c(0,1)))]
# inds<-tmp=="[Term]"
# runs<-rle(inds)
# runs$values<-rep(1:(length(runs$values)/2),each=2)
# inds<-inverse.rle(runs)
# tmp<-split(tmp,inds)
# ids<-lapply(tmp,function(ii) gsub("^id\\: ","",ii[grepl("^id\\:",ii)]))
# any(lengths(ids)>1) #no double ids
# pas<-lapply(tmp,
#             function(ii){
#               tmp<-gsub("^is_a\\: ","",ii[grepl("^is_a\\:",ii)])
#               tmp<-substr(tmp,1,regexpr("\\!",tmp)-2)
#             })
# lookup<-setNames(pas,ids)
# 
# #parent- getting
# #perhaps slower than necessary, but straight-forward and reasonably fast
# all.parents<-setNames(apply(do.call(cbind,lapply(allgo[,-1],strsplit,split="\\|")),1,unlist,use.names=FALSE),allgo[,1])
# for(i in seq_along(all.parents)){
#   qq<-all.parents[[i]]
#   while(length(qq)>1){
#     qq<-unlist(lookup[match(qq,names(lookup))],use.names=FALSE)
#     all.parents[[i]]<-unique(c(all.parents[[i]],qq))
#   }
#   if(!(i%%1000)) cat(i,"\n")
# }
# #just for comparison to make sure it looks right...
# og<-setNames(apply(do.call(cbind,lapply(allgo[,-1],strsplit,split="\\|")),1,unlist,use.names=FALSE),allgo[,1])
# plot(lengths(all.parents)~lengths(og));abline(0,1)
# any(lengths(all.parents)<lengths(og))
# #welp, the number of parents is always equal to or greater than the number of children, so that's promising!
# #some empty list elements, but they seem to correspond to rows without go terms in "allgo":
# all((!lengths(og))==(!nchar(allgo[,2])&!nchar(allgo[,3])))
# 
# all(allgo[,1]==names(all.parents))
# allgo$all.pars<-unlist(lapply(all.parents,paste,collapse="|"),use.names=FALSE)
# 
# write.csv(allgo, "Vriparia-functional-annotations-withparents.tsv", sep='\t')

rn1 <- paste(allgo[,1], sep="")
#gene2GO <- paste(allgo$Arabidopsis_GO_terms, allgo$PFAM_GO_terms, sep="|")
#gene2GO = allgo[,-1]
gene2GO <- paste(allgo$all.pars)
names(gene2GO)<-rn1
gene2GO<-strsplit(gene2GO,"\\|")

# remove any candidate genes without GO annotation
candidate_list <- degenes$Gene
keep = candidate_list %in% allgo$gene
keep =which(keep==TRUE)
candidate_list=candidate_list[keep]

# make named factor showing which genes are of interest
bg_genes <- allgo$gene
geneList=factor(as.integer(bg_genes %in% candidate_list))
names(geneList)= bg_genes


## Set up as topGOdata with BP ontology
GODATA <- new("topGOdata",
              description = "conatia vs Leaf in 588con",
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
all_res_BP_con=GenTable(GODATA, weightFisher=resultFisher.BP, orderBy='score', topNodes=length(allGO_BP))
sign_res_BP_con <- all_res_BP_con[all_res_BP_con$weightFisher<0.05,]

#performing BH correction on our p values
p.adj=round(p.adjust(all_res_BP_con$weightFisher,method="BH"),digits = 4)

# create the file with all the statistics from GO analysis
all_res_BP_con_final=cbind(all_res_BP_con,p.adj)
all_res_BP_con_final=all_res_BP_con_final[order(all_res_BP_con_final$p.adj),]

sign_res_BP_con <- all_res_BP_con[all_res_BP_con$weightFisher<0.05,]
sign_res_BP_con_final <- all_res_BP_con_final[all_res_BP_con_final$p.adj<0.05,]

## Set up as topGOdata with MF ontology
GODATA <- new("topGOdata",
              description = "conatia vs Leaf in 588con",
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
all_res_MF_con=GenTable(GODATA, weightFisher=resultFisher.MF, orderBy='score', topNodes=length(allGO_MF))
sign_res_MF_con <- all_res_MF_con[all_res_MF_con$weightFisher<0.05,]

#performing BH correction on our p values
p.adj=round(p.adjust(all_res_MF_con$weightFisher,method="BH"),digits = 4)

# create the file with all the statistics from GO analysis
all_res_MF_con_final=cbind(all_res_MF_con,p.adj)
all_res_MF_con_final=all_res_MF_con_final[order(all_res_MF_con_final$p.adj),]

sign_res_MF_con <- all_res_MF_con[all_res_MF_con$weightFisher<0.05,]
sign_res_MF_con_final <- all_res_MF_con_final[all_res_MF_con_final$p.adj<0.05,]

## Set up as topGOdata with CC ontology
GODATA <- new("topGOdata",
              description = "conatia vs Leaf in 588con",
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
all_res_CC_con=GenTable(GODATA, weightFisher=resultFisher.CC, orderBy='score', topNodes=length(allGO_CC))
sign_res_CC_con <- all_res_CC_con[all_res_CC_con$weightFisher<0.05,]


#performing BH correction on our p values
p.adj=round(p.adjust(all_res_CC_con$weightFisher,method="BH"),digits = 4)

# create the file with all the statistics from GO analysis
all_res_CC_con_final=cbind(all_res_CC_con,p.adj)
all_res_CC_con_final=all_res_CC_con_final[order(all_res_CC_con_final$p.adj),]

sign_res_CC_con <- all_res_CC_con[all_res_CC_con$weightFisher<0.05,]
sign_res_CC_con_final <- all_res_CC_con_final[all_res_CC_con_final$p.adj<0.05,]
