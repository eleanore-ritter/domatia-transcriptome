## NOTE: This code will need to be gone through extensively to make sure we keep only relevant analyses

# Set working directory and load necessary packages
setwd("C:/Users/rittere5/OneDrive - Michigan State University/Vitis-domatia/")
#setwd("C:/Users/elean/OneDrive - Michigan State University/Vitis-domatia/")

library(dplyr)
library(tidyr)
library(ggplot2)

# Read in differentially expressed genes
dom710 <- read.csv("differentially_expressed_genes/DE_genes_Domatia_V_Leaf_588710.csv")
dom711 <- read.csv("differentially_expressed_genes/DE_genes_Domatia_V_Leaf_588711.csv")

dom710.upreg <- dom710[dom710$log2FoldChange>0 ,]
dom710.downreg <- dom710[dom710$log2FoldChange<0 ,]

dom711.upreg <- dom711[dom711$log2FoldChange>0 ,]
dom711.downreg <- dom711[dom711$log2FoldChange<0 ,]

# Read in EFN transport genes
efn <- read.csv("efn-genes-all.csv") # Looking at all genes
#efn <- read.csv("efn-genes-transporters.csv") # Looking at transporters specifically
efn$A.thaliana.ID <- gsub("\\.[0-9]", "", efn$A.thaliana.ID)

# Filter out genes differentially expressed in foliar nectaries compared to control
foliar.degenes <- efn[efn$DE..FO.post.vs.FO.post.C=="Y" | efn$DE..FO.pre.vs.FO.pre.C=="Y" | efn$DE..FO.sec.vs.FO.sec.C=="Y" ,]
foliar.degenes$A.thaliana.ID <- gsub("\\.[0-9]", "", foliar.degenes$A.thaliana.ID)
foliar.upreg <- foliar.degenes[foliar.degenes$DE..FO.post.vs.FO.post.C.Up.or.Down=="Up" & foliar.degenes$DE..FO.pre.vs.FO.pre.C.Up.or.Down=="Up" &
                 foliar.degenes$DE..FO.sec.vs.FO.sec.C.Up.or.Down=="Up" ,]

foliar.downreg <- foliar.degenes[foliar.degenes$DE..FO.post.vs.FO.post.C.Up.or.Down=="Down" & foliar.degenes$DE..FO.pre.vs.FO.pre.C.Up.or.Down=="Down" &
                                 foliar.degenes$DE..FO.sec.vs.FO.sec.C.Up.or.Down=="Down" ,]

foliarsec.upreg <- foliar.degenes[foliar.degenes$DE..FO.sec.vs.FO.sec.C.Up.or.Down=="Up" ,]
foliarpre.upreg <- foliar.degenes[foliar.degenes$DE..FO.pre.vs.FO.pre.C.Up.or.Down=="Up" ,]

foliar.anystage.upreg <- foliar.degenes[foliar.degenes$DE..FO.post.vs.FO.post.C.Up.or.Down=="Up" | foliar.degenes$DE..FO.pre.vs.FO.pre.C.Up.or.Down=="Up" |
                                 foliar.degenes$DE..FO.sec.vs.FO.sec.C.Up.or.Down=="Up" ,]  

# Get overlap for each genotype individually
shared710.FOupreg <- merge(foliar.upreg, dom710.upreg, by.x="A.thaliana.ID",by.y="ensembl_gene_id")
shared711.FOlupreg <- merge(foliar.upreg, dom711.upreg, by.x="A.thaliana.ID",by.y="ensembl_gene_id")

sharedsec710.FOlupreg <- merge(foliarsec.upreg, dom710.upreg, by.x="A.thaliana.ID",by.y="ensembl_gene_id")
sharedpre710.FOupreg <- merge(foliarpre.upreg, dom710.upreg, by.x="A.thaliana.ID",by.y="ensembl_gene_id")

shared710.anystage.FOupreg <-  merge(foliar.anystage.upreg, dom710.upreg, by.x="A.thaliana.ID",by.y="ensembl_gene_id")
shared711.anystage.FOupreg <-  merge(foliar.anystage.upreg, dom711.upreg, by.x="A.thaliana.ID",by.y="ensembl_gene_id")

# Overlap for both genotypes shared
dom.bothgeno.upreg <- dom710.upreg[dom710.upreg$Gene %in% dom711.upreg$Gene,]

sharedboth.anystage.upreg <- merge(foliar.anystage.upreg, dom.bothgeno.upreg, by.x="A.thaliana.ID",by.y="ensembl_gene_id")

# Repeat for ANY foliar, bracteal, and circumbracteal upregulated genes at ANY stage
efn.degenes <- efn[efn$DE..FO.post.vs.FO.post.C=="Y" | efn$DE..FO.pre.vs.FO.pre.C=="Y" | efn$DE..FO.sec.vs.FO.sec.C=="Y" | 
                     efn$DE..B.post.vs.B.post.C=="Y" | efn$DE..B.pre.vs.B.pre.C=="Y" | efn$DE..B.sec.vs.B.sec.C=="Y" |
                     efn$DE..C.post.vs.C.post.C=="Y" | efn$DE..C.pre.vs.C.pre.C=="Y" | efn$DE..C.sec.vs.C.sec.C=="Y" ,]
efn.degenes.upreg <- efn.degenes[efn.degenes$DE..FO.post.vs.FO.post.C.Up.or.Down=="Up" | efn.degenes$DE..FO.pre.vs.FO.pre.C.Up.or.Down=="Up" |
                                          efn.degenes$DE..FO.sec.vs.FO.sec.C.Up.or.Down=="Up" |
                                          efn.degenes$DE..B.post.vs.B.post.C.Up.or.Down=="Up" | efn.degenes$DE..B.pre.vs.B.pre.C.Up.or.Down=="Up" |
                                          efn.degenes$DE..B.sec.vs.B.sec.C.Up.or.Down=="Up" |
                                          efn.degenes$DE..C.post.vs.C.post.C.Up.or.Down=="Up" | efn.degenes$DE..C.pre.vs.C.pre.C.Up.or.Down=="Up" |
                                          efn.degenes$DE..C.sec.vs.C.sec.C.Up.or.Down=="Up" ,]  
shared710.anystage.ALLupreg <-  merge(efn.degenes.upreg, dom710.upreg, by.x="A.thaliana.ID",by.y="ensembl_gene_id")
shared711.anystage.ALLupreg <-  merge(efn.degenes.upreg, dom711.upreg, by.x="A.thaliana.ID",by.y="ensembl_gene_id")
sharedboth.anystage.ALLupreg <- merge(efn.degenes.upreg, dom.bothgeno.upreg, by.x="A.thaliana.ID",by.y="ensembl_gene_id")

test <-sharedboth.anystage.ALLupreg[sharedboth.anystage.ALLupreg$A.thaliana.ID=="AT1G22710",]
