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
col.line <- c("#1158f6") 

lline1 <- c(2.6 ,2.6)
lline2 <- c(0, 0.125)
lline <- data.frame(lline1, lline2)
colnames(lline) <- c("Leaf_Width", "Radius")
lline$Genotype <- c("line")

uline1 <- c(5.8 ,5.8)
uline2 <- c(0, 0.125)
uline <- data.frame(uline1, uline2)
colnames(uline) <- c("Leaf_Width", "Radius")
uline$Genotype <- c("line")

ggplot(data=data, aes(x=Leaf_Width, y=Radius, color=Genotype)) +
  geom_point(size = 2, alpha = 0.9) +
  scale_color_manual(values = c(col.710, col.711, col.line),
                     labels = c("588710 (SDG)", "588711 (LDG)", NULL)) +
  theme_classic() +
  geom_smooth(method='lm') +
  scale_y_continuous(limits=c(0,0.125), expand = c(0,0) ) +
  scale_x_continuous(expand = c(0,0) ) +
  geom_path(data=lline, linewidth=2) +
  geom_path(data=uline, linewidth=2) +
  xlab("Leaf Width (cm)") +
  ylab("Domatium Radius (mm)") +
  theme(axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12)) 

########################## STATS ##########################
lmtest<-lm(Leaf_Width~Radius, data=data)
summary(test)

midage <- data[data$Leaf_Width>2.6 & data$Leaf_Width<5.8 ,1]
oldage <- data[data$Leaf_Width>5.8 ,1]
ttest <- t.test(midage, oldage)
