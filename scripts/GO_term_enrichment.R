# Set working directory and load necessary packages
setwd("C:/Users/rittere5/OneDrive - Michigan State University/Vitis-domatia/")
#setwd("C:/Users/elean/OneDrive - Michigan State University/Vitis-domatia/")

library("topGO")

# Load data
all <- read.csv("Vriparia-functional-annotations-V2.tsv", sep='\t')
degenes <- read.csv("DE_genes_Domatia_V_Leaf_588710.csv")

# Modify data

## Modify GO term file
allgo <- all[,c(2,5)]
rn1 <- paste(allgo[,1], sep="")
gene2GO = allgo[,-1]
names(gene2GO)<-rn1
gene2GO<-strsplit(gene2GO,"\\|")

## Modify DE gene list
matches<-setNames(match(degenes$Gene,allgo$gene,nomatch=0),degenes$Gene)
matches[matches>0]<-1

# Set up GO data
## Set up function
ft <- function(x){
  return(x>0)
}

## Set up as topGOdata
GODATA <- new("topGOdata",
              description = "Domatia vs Leaf in 588710",
              ontology = "CC",
              allGenes = matches,
              gene2GO = gene2GO,
              nodeSize = 10,
              geneSel = ft,
              annotationFun = annFUN.gene2GO
              )

# Run GO term enrichment analysis

## Run with Fisher's exact test
resultFisher <- runTest(GODATA, algorithm = "classic", statistic = "fisher")
resultFisher

## Run with Kolmogorov-Smirnov test
resultKS <- runTest(GODATA, algorithm = "classic", statistic = "ks") # Classic method
resultKS

resultKS.elim <- runTest(GODATA, algorithm = "elim", statistic = "ks") # Elimination method
resultKS.elim
