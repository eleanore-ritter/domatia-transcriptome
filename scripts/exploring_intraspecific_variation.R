# Set working directory and load necessary packages
#setwd("C:/Users/rittere5/OneDrive - Michigan State University/Vitis-domatia/")
setwd("C:/Users/elean/OneDrive - Michigan State University/Vitis-domatia/")

library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(ggcorrplot)
library(viridis)

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
# to highly expressed in domatia tissue, and alternate genotype just has leaf expression=domatia expression
