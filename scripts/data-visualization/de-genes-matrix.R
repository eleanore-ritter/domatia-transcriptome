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

############################# WITHOUT genes unique to genotypes #############################
# Make matrix of overlapping genes - WITHOUT genes unique to each genotype

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

############################# WITH genes unique to genotypes #############################
# Make matrix of overlapping genes - WITH genes unique to each genotype

mat = matrix(
  c(1447, NA, NA, NA, NA, NA, 909, 909, NA, NA, NA, NA, 538, 0, 759, NA, NA, NA, 0, 0, 221, 221, NA, NA,
    13, 9, 9, 5, 19, NA, 459, 319, 227, 87, 16, 2888),
  nrow = 6,  
  ncol = 6,        
  byrow = TRUE         
)


rownames(mat) = c("C vs D - 588710", "C vs D - 588710 Unique", "C vs D - 588711", "C vs D - 588711 Unique", "D vs D", "C vs C")

colnames(mat) = c("C vs D - 588710", "C vs D - 588710 Unique", "C vs D - 588711", "C vs D - 588711 Unique", "D vs D", "C vs C")

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

