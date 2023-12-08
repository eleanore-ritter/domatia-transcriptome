# Set working directory and load necessary packages
setwd("C:/Users/rittere5/OneDrive - Michigan State University/Vitis-domatia/")
#setwd("C:/Users/elean/OneDrive - Michigan State University/Vitis-domatia/")

library(ggplot2)
library(cowplot)
library(ggsignif)
library(tidyr)

########################## SET UP DATA IN R ##########################
# Load raw data
raw.data <- read.csv("domatium ontogeny data.csv")

# Extract out our two genotypes
gen1 <- raw.data[grep("588710",raw.data$plant_code),]
gen2 <- raw.data[grep("588711",raw.data$plant_code),]

# Merge all radius and density data
gen1.temp1 <- data.frame(gen1$radius_1, gen1$radius_2, gen1$radius_3, gen1$radius_4)
gen1.temp2 <- as.data.frame(unlist(gen1.temp1))
gen1.temp3 <- as.data.frame(unlist(data.frame(gen1$density_1, gen1$density_2, gen1$density_3, gen1$density_4)))
gen1.temp3$genotype <- c("588710")
gen1.temp4 <- as.data.frame(unlist(data.frame(gen1$leaf)))
myLetters <- LETTERS[1:26]
gen1.temp5 <- as.factor(match(gen1.temp4$`unlist(data.frame(gen1$leaf))`, myLetters))
gen1.temp6 <- as.data.frame(unlist(data.frame(gen1$leaf_width)))
gen1.fin <- cbind(gen1.temp2, gen1.temp3, gen1.temp5, gen1.temp6)
colnames(gen1.fin) <- c("Radius", "Density", "Genotype", "Leaf", "Leaf_Width")

gen2.temp1 <- data.frame(gen2$radius_1, gen2$radius_2, gen2$radius_3, gen2$radius_4)
gen2.temp2 <- as.data.frame(unlist(gen2.temp1))
gen2.temp3 <- as.data.frame(unlist(data.frame(gen2$density_1, gen2$density_2, gen2$density_3, gen2$density_4)))
gen2.temp3$genotype <- c("588711")
gen2.temp4 <- as.data.frame(unlist(data.frame(gen2$leaf)))
myLetters <- LETTERS[1:26]
gen2.temp5 <- as.factor(match(gen2.temp4$`unlist(data.frame(gen2$leaf))`, myLetters))
gen2.temp6 <- as.data.frame(unlist(data.frame(gen2$leaf_width)))
gen2.fin <- cbind(gen2.temp2, gen2.temp3, gen2.temp5, gen2.temp6)
colnames(gen2.fin) <- c("Radius", "Density", "Genotype", "Leaf", "Leaf_Width")

data <- rbind(gen1.fin, gen2.fin)
#data$Genotype <- gsub("588710", "Little", data$Genotype)
#data$Genotype <- gsub("588711", "Big", data$Genotype)
data$Genotype <- as.factor(data$Genotype)
data$Radius <- as.numeric(data$Radius)
data$Density <- as.numeric(data$Density)
data$Leaf_Width <- as.numeric(data$Leaf_Width)

########################## PLOT PHENOTYPE DATA ##########################

col.711 <- c("#024F4A")
col.710 <- c("#05B384")

ggplot(data, aes(x=Leaf_Width, y=Radius, color=Genotype)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(values = c(col.710, col.711)) +
  theme_classic() +
  geom_smooth(method='lm')

########################## STATS ##########################

