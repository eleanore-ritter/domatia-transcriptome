# Set working directory and load necessary packages
setwd("C:/Users/rittere5/OneDrive - Michigan State University/Vitis-domatia/")
#setwd("C:/Users/elean/OneDrive - Michigan State University/Vitis-domatia/")

library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(ggcorrplot)
library(viridis)
library("ggVennDiagram")

# Load data
concon <- read.csv("differentially_expressed_genes/DE_genes_Control_710_v_711.csv")
row.names(concon) <- concon$Gene

domdom <- read.csv("differentially_expressed_genes/DE_genes_Domatia_710_v_711.csv")
row.names(domdom) <- domdom$Gene

domcon.710 <- read.csv("differentially_expressed_genes/DE_genes_Domatia_V_Leaf_588710.csv")
row.names(domcon.710) <- domcon.710$Gene

domcon.711 <- read.csv("differentially_expressed_genes/DE_genes_Domatia_V_Leaf_588711.csv")
row.names(domcon.711) <- domcon.711$Gene

############################# Plot venn diagram #############################
x <- list(
  A = rownames(domcon.710),
  B = rownames(concon),
  C = rownames(domdom),
  D = rownames(domcon.711)
)

venn <- Venn(x)
data <- process_data(venn)
ggplot() +
  # 1. region count layer
  geom_sf(aes(fill = count), data = venn_region(data)) +
  # 2. set edge layer
  geom_sf(color = "black", data = venn_setedge(data), show.legend = FALSE) +
  # 3. set label layer
  geom_sf_text(aes(label = c("588710:\nC vs D", "Control vs Control", "Domatia vs Domatia", "588711:\nC vs D")),data = venn_setlabel(data), fontface = "bold", size = 4.8) +
  # 4. region label layer
  geom_sf_label(aes(label = count), data = venn_region(data), alpha = 0, label.size = 0) +
  theme_void() +
  scale_fill_gradient2(limit = c(0,2500), low = "#fcfdbf", high =  "#ca3e72", mid = "#fe9f6d", midpoint = 1250) +
  coord_sf(default_crs = NULL) +
  labs(fill = "Differentially\nexpressed genes\n")+
  guides(fill = guide_colourbar(title.position = "top", frame.colour = "black", ticks.colour = "black", draw.ulim = FALSE, draw.llim = FALSE))