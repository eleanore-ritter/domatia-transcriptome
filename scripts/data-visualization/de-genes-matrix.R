# Set working directory and load necessary packages
setwd("C:/Users/rittere5/OneDrive - Michigan State University/Vitis-domatia/")
#setwd("C:/Users/elean/OneDrive - Michigan State University/Vitis-domatia/")

library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(ggcorrplot)
library(viridis)

# Load data
concon <- read.csv("differentially_expressed_genes/DE_genes_Control_710_v_711.csv")
row.names(concon) <- concon$Gene

domdom <- read.csv("differentially_expressed_genes/DE_genes_Domatia_710_v_711.csv")
row.names(domdom) <- domdom$Gene

domcon.710 <- read.csv("differentially_expressed_genes/DE_genes_Domatia_V_Leaf_588710.csv")
row.names(domcon.710) <- domcon.710$Gene

domcon.711 <- read.csv("differentially_expressed_genes/DE_genes_Domatia_V_Leaf_588711.csv")
row.names(domcon.711) <- domcon.711$Gene



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



# Make matrix of overlapping gene

mat = matrix(
  c(1447, 538, 13, 459, 538, 759, 9, 227, 13, 9, 19, 16, 459, 227, 16, 2888),
  nrow = 4,  
  ncol = 4,        
  byrow = TRUE         
)

mat = matrix(
  c(1447, NA, NA, NA, 538, 759, NA, NA, 13, 9, 19, NA, 459, 227, 16, 2888),
  nrow = 4,  
  ncol = 4,        
  byrow = TRUE         
)


rownames(mat) = c("C vs D - 588710", "C vs D - 588711", "D vs D", "C vs C")

colnames(mat) = c("C vs D - 588710", "C vs D - 588711", "D vs D", "C vs C")

print(mat)



# Plot matrix

p <- ggcorrplot(mat, lab = TRUE, outline.col = "black", tl.srt = 45)
p + 
  scale_fill_gradient2(limit = c(0,3000), low = "#fcfdbf", high =  "#ca3e72", mid = "#fe9f6d", midpoint = 1500) + 
  theme(axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"),
        legend.position = "top",
        legend.key.width = unit(1, "cm")) +
  labs(fill = "Number of DE Genes") +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5, frame.colour = "black", ticks.colour = "black", draw.ulim = FALSE, draw.llim = FALSE))

