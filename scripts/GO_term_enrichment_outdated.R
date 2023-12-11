# Set working directory and load necessary packages
#setwd("C:/Users/rittere5/OneDrive - Michigan State University/Vitis-domatia/")
setwd("C:/Users/elean/OneDrive - Michigan State University/Vitis-domatia/")

library("topGO")
library(dplyr)
library(tidyr)
library(GOSim)
library(GO.db)

######################## GO TERM ENRICHMENT WITH DOMATIA V LEAF IN 588710 ######################## 

# Load data

## Load functional annotations
all <- read.csv("Vriparia-functional-annotations-TAIR10-mod.tsv", sep='\t')

## Add gene data to all if needed
dict <- read.csv("cds-to-gene.tsv", sep = '\t')
colnames(dict) <- c("Transcript", "gene")
dict$Transcript <- gsub("\\.[0-9]", "", dict$Transcript)
all <- merge(all, dict, by="Transcript", all = TRUE)
all <- all[,c(1, 8, 3:7)]
all <- all %>% drop_na(gene)
all <- all[!duplicated(all$gene),] # Get rid of duplicates - all rows are identical (transcripts are sometimes different, that is all)

## Proving that duplicated rows are identical for all gene names and GO terms as a sanity check

test <- all
check<-do.call(rbind,tapply(1:nrow(test),
                            test$gene,
                            function(ii) lengths(lapply(test,function(jj) unique(jj[ii])))))
apply(check,2,function(ii) any(ii>2))

## Load differentially expressed genes
degenes <- read.csv("differentially_expressed_genes/DE_genes_Domatia_V_Leaf_588710.csv")


# Modify data

## Modify GO term file
allgo <- all[,c(2,4,6)]
foo<-function(x){
  inds<-allgo$gene==x
  tmp1<-unique(unlist(strsplit(allgo$Arabidopsis_GO_terms[inds],"\\|"),use.names=FALSE))
  tmp1<-tmp1[nzchar(tmp1)&!is.na(tmp1)]
  tmp2<-unique(unlist(strsplit(allgo$PFAM_GO_terms[inds],"\\|"),use.names=FALSE))
  tmp2<-tmp2[nzchar(tmp2)&!is.na(tmp2)]
  data.frame("gene"=x,
             "Arabidopsis_GO_terms"=paste(tmp1,collapse="|"),
             "PFAM_GO_terms"=paste(tmp2,collapse="|"))
}
allgo <- do.call(rbind,
                 lapply(unique(allgo$gene),foo))


## Add back in parent GO terms
gotest <- c("GO:0000812", "GO:0000742")

xx<- as.list(GOBPPARENTS)
if(length(gotest) > 0){
  # Get the parent GO IDs for the first elents of xx
  goids <- gotest[[1]]
  # Find out the GO terms for the first parent goid
  GOID(GOTERM[[goids[1]]])
  Term(GOTERM[[goids[1]]])
  Synonym(GOTERM[[goids[1]]])
  Secondary(GOTERM[[goids[1]]])
  Definition(GOTERM[[goids[1]]])
  Ontology(GOTERM[[goids[1]]])
}
GOBPPARENTS(gotest)

rn1 <- paste(allgo[,1], sep="")
gene2GO <- paste(allgo$Arabidopsis_GO_terms, allgo$PFAM_GO_terms, sep="|")
#gene2GO = allgo[,-1]
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

## Set up as topGOdata with BP ontology
GODATA <- new("topGOdata",
              description = "Domatia vs Leaf in 588710",
              ontology = "BP",
              allGenes = matches,
              gene2GO = gene2GO,
              geneSel = ft,
              annotationFun = annFUN.gene2GO,
              )

# Run GO term enrichment analysis

## Run with Fisher's exact test
resultFisher.BP <- runTest(GODATA, algorithm = "classic", statistic = "fisher")
resultFisher.BP

## Run with Kolmogorov-Smirnov test
resultKS.BP710 <- runTest(GODATA, algorithm = "classic", statistic = "ks") # Classic method
resultKS.BP710
pvalKS.BP710 <- as.data.frame(score(resultKS.BP710))

# resultKS.elim.BP <- runTest(GODATA, algorithm = "elim", statistic = "ks") # Elimination method
# resultKS.elim.BP

## Set up as topGOdata with MF ontology
GODATA <- new("topGOdata",
              description = "Domatia vs Leaf in 588710",
              ontology = "MF",
              allGenes = matches,
              gene2GO = gene2GO,
              geneSel = ft,
              annotationFun = annFUN.gene2GO
)

# Run GO term enrichment analysis

# ## Run with Fisher's exact test
# resultFisher.MF <- runTest(GODATA, algorithm = "classic", statistic = "fisher")
# resultFisher.MF

## Run with Kolmogorov-Smirnov test
resultKS.MF710 <- runTest(GODATA, algorithm = "classic", statistic = "ks") # Classic method
resultKS.MF710
pvalKS.MF710 <- as.data.frame(score(resultKS.MF710))

# resultKS.elim.MF <- runTest(GODATA, algorithm = "elim", statistic = "ks") # Elimination method
# resultKS.elim.MF

## Set up as topGOdata with CC ontology
GODATA <- new("topGOdata",
              description = "Domatia vs Leaf in 588710",
              ontology = "CC",
              allGenes = matches,
              gene2GO = gene2GO, 
              geneSel = ft,
              annotationFun = annFUN.gene2GO
)

# Run GO term enrichment analysis

## Run with Fisher's exact test
# resultFisher.CC <- runTest(GODATA, algorithm = "classic", statistic = "fisher")
# resultFisher.CC

## Run with Kolmogorov-Smirnov test
resultKS.CC710 <- runTest(GODATA, algorithm = "classic", statistic = "ks") # Classic method
resultKS.CC710
pvalKS.CC710 <- as.data.frame(score(resultKS.CC710))

# resultKS.elim.CC <- runTest(GODATA, algorithm = "elim", statistic = "ks") # Elimination method

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
degenes <- read.csv("DE_genes_Domatia_V_Leaf_588711.csv")

# Modify data

## Modify GO term file
allgo <- all[,c(2,7)]
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

## Set up as topGOdata with BP ontology
GODATA <- new("topGOdata",
              description = "Domatia vs Leaf in 588711",
              ontology = "BP",
              allGenes = matches,
              gene2GO = gene2GO,
              geneSel = ft,
              annotationFun = annFUN.gene2GO
)

# Run GO term enrichment analysis

## Run with Fisher's exact test
# resultFisher.BP <- runTest(GODATA, algorithm = "classic", statistic = "fisher")
# resultFisher.BP

## Run with Kolmogorov-Smirnov test
resultKS.BP711 <- runTest(GODATA, algorithm = "classic", statistic = "ks") # Classic method
resultKS.BP711
pvalKS.BP711 <- as.data.frame(score(resultKS.BP711))

# resultKS.elim.BP <- runTest(GODATA, algorithm = "elim", statistic = "ks") # Elimination method
# resultKS.elim.BP

## Set up as topGOdata with MF ontology
GODATA <- new("topGOdata",
              description = "Domatia vs Leaf in 588710",
              ontology = "MF",
              allGenes = matches,
              gene2GO = gene2GO,
              geneSel = ft,
              annotationFun = annFUN.gene2GO
)

# Run GO term enrichment analysis

# ## Run with Fisher's exact test
# resultFisher.MF <- runTest(GODATA, algorithm = "classic", statistic = "fisher")
# resultFisher.MF

## Run with Kolmogorov-Smirnov test
resultKS.MF711 <- runTest(GODATA, algorithm = "classic", statistic = "ks") # Classic method
resultKS.MF711
pvalKS.MF711 <- as.data.frame(score(resultKS.MF711))

# resultKS.elim.MF <- runTest(GODATA, algorithm = "elim", statistic = "ks") # Elimination method
# resultKS.elim.MF

## Set up as topGOdata with CC ontology
GODATA <- new("topGOdata",
              description = "Domatia vs Leaf in 588710",
              ontology = "CC",
              allGenes = matches,
              gene2GO = gene2GO, 
              geneSel = ft,
              annotationFun = annFUN.gene2GO
)

# Run GO term enrichment analysis

## Run with Fisher's exact test
# resultFisher.CC <- runTest(GODATA, algorithm = "classic", statistic = "fisher")
# resultFisher.CC

## Run with Kolmogorov-Smirnov test
resultKS.CC711 <- runTest(GODATA, algorithm = "classic", statistic = "ks") # Classic method
resultKS.CC711
pvalKS.CC711 <- as.data.frame(score(resultKS.CC711))

# resultKS.elim.CC <- runTest(GODATA, algorithm = "elim", statistic = "ks") # Elimination method
# resultKS.elim.CC
# resultKS.elim.CC

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
degenes <- read.csv("DE_genes_Domatia_710_v_711.csv")

# Modify data

## Modify GO term file
allgo <- all[,c(2,7)]
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

## Set up as topGOdata with BP ontology
GODATA <- new("topGOdata",
              description = "Domatia vs Domatia",
              ontology = "BP",
              allGenes = matches,
              gene2GO = gene2GO,
              geneSel = ft,
              annotationFun = annFUN.gene2GO
)

# Run GO term enrichment analysis

## Run with Fisher's exact test
# resultFisher.BP <- runTest(GODATA, algorithm = "classic", statistic = "fisher")
# resultFisher.BP

## Run with Kolmogorov-Smirnov test
resultKS.BPDOM <- runTest(GODATA, algorithm = "classic", statistic = "ks") # Classic method
resultKS.BPDOM
pvalKS.BPDOM <- as.data.frame(score(resultKS.BPDOM))

# resultKS.elim.BP <- runTest(GODATA, algorithm = "elim", statistic = "ks") # Elimination method
# resultKS.elim.BP

## Set up as topGOdata with MF ontology
GODATA <- new("topGOdata",
              description = "Domatia vs Domatia",
              ontology = "MF",
              allGenes = matches,
              gene2GO = gene2GO,
              geneSel = ft,
              annotationFun = annFUN.gene2GO
)

# Run GO term enrichment analysis

# ## Run with Fisher's exact test
# resultFisher.MF <- runTest(GODATA, algorithm = "classic", statistic = "fisher")
# resultFisher.MF

## Run with Kolmogorov-Smirnov test
resultKS.MFDOM <- runTest(GODATA, algorithm = "classic", statistic = "ks") # Classic method
resultKS.MFDOM
pvalKS.MFDOM <- as.data.frame(score(resultKS.MFDOM))

# resultKS.elim.MF <- runTest(GODATA, algorithm = "elim", statistic = "ks") # Elimination method
# resultKS.elim.MF

## Set up as topGOdata with CC ontology
GODATA <- new("topGOdata",
              description = "Domatia vs Domatia",
              ontology = "CC",
              allGenes = matches,
              gene2GO = gene2GO, 
              geneSel = ft,
              annotationFun = annFUN.gene2GO
)

# Run GO term enrichment analysis

## Run with Fisher's exact test
# resultFisher.CC <- runTest(GODATA, algorithm = "classic", statistic = "fisher")
# resultFisher.CC

## Run with Kolmogorov-Smirnov test
resultKS.CCDOM <- runTest(GODATA, algorithm = "classic", statistic = "ks") # Classic method
resultKS.CCDOM
pvalKS.CCDOM <- as.data.frame(score(resultKS.CCDOM))

# resultKS.elim.CC <- runTest(GODATA, algorithm = "elim", statistic = "ks") # Elimination method
# resultKS.elim.CC
# resultKS.elim.CC

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
degenes <- read.csv("DE_genes_Control_710_v_711.csv")

# Modify data

## Modify GO term file
allgo <- all[,c(2,7)]
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

## Set up as topGOdata with BP ontology
GODATA <- new("topGOdata",
              description = "Control v Control",
              ontology = "BP",
              allGenes = matches,
              gene2GO = gene2GO,
              geneSel = ft,
              annotationFun = annFUN.gene2GO
)

# Run GO term enrichment analysis

## Run with Fisher's exact test
# resultFisher.BP <- runTest(GODATA, algorithm = "classic", statistic = "fisher")
# resultFisher.BP

## Run with Kolmogorov-Smirnov test
resultKS.BPCONTROL <- runTest(GODATA, algorithm = "classic", statistic = "ks") # Classic method
resultKS.BPCONTROL
pvalKS.BPCONTROL <- as.data.frame(score(resultKS.BPCONTROL))

# resultKS.elim.BP <- runTest(GODATA, algorithm = "elim", statistic = "ks") # Elimination method
# resultKS.elim.BP

## Set up as topGOdata with MF ontology
GODATA <- new("topGOdata",
              description = "Control v Control",
              ontology = "MF",
              allGenes = matches,
              gene2GO = gene2GO,
              geneSel = ft,
              annotationFun = annFUN.gene2GO
)

# Run GO term enrichment analysis

# ## Run with Fisher's exact test
# resultFisher.MF <- runTest(GODATA, algorithm = "classic", statistic = "fisher")
# resultFisher.MF

## Run with Kolmogorov-Smirnov test
resultKS.MFCONTROL <- runTest(GODATA, algorithm = "classic", statistic = "ks") # Classic method
resultKS.MFCONTROL
pvalKS.MFCONTROL <- as.data.frame(score(resultKS.MFCONTROL))

# resultKS.elim.MF <- runTest(GODATA, algorithm = "elim", statistic = "ks") # Elimination method
# resultKS.elim.MF

## Set up as topGOdata with CC ontology
GODATA <- new("topGOdata",
              description = "Control v Control",
              ontology = "CC",
              allGenes = matches,
              gene2GO = gene2GO, 
              geneSel = ft,
              annotationFun = annFUN.gene2GO
)

# Run GO term enrichment analysis

## Run with Fisher's exact test
# resultFisher.CC <- runTest(GODATA, algorithm = "classic", statistic = "fisher")
# resultFisher.CC

## Run with Kolmogorov-Smirnov test
resultKS.CCCONTROL <- runTest(GODATA, algorithm = "classic", statistic = "ks") # Classic method
resultKS.CCCONTROL
pvalKS.CCCONTROL <- as.data.frame(score(resultKS.CCCONTROL))

# resultKS.elim.CC <- runTest(GODATA, algorithm = "elim", statistic = "ks") # Elimination method
# resultKS.elim.CC
# resultKS.elim.CC