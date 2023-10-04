# Set working directory and load necessary packages
setwd("C:/Users/rittere5/OneDrive - Michigan State University/Vitis-domatia/")
#setwd("C:/Users/elean/OneDrive - Michigan State University/Vitis-domatia/")

library(ggplot2)
library(cowplot)
library(ggsignif)
library(tidyr)

########################## SET UP DATA IN R ##########################
# Load raw data
raw.data <- read.csv("domatia_trichomes_2022_final.csv")

# Extract out our two genotypes
gen1 <- raw.data[grep("588710",raw.data$plant_code),]
gen2 <- raw.data[grep("588711",raw.data$plant_code),]

# Merge all radius and density data
gen1.temp1 <- data.frame(gen1$radius_1, gen1$radius_2, gen1$radius_3, gen1$radius_4)
gen1.temp2 <- as.data.frame(unlist(gen1.temp1))
gen1.temp3 <- as.data.frame(unlist(data.frame(gen1$density_1, gen1$density_2, gen1$density_3, gen1$density_4)))
gen1.temp3$genotype <- c("588710")
gen1.fin <- cbind(gen1.temp2, gen1.temp3)
colnames(gen1.fin) <- c("Radius", "Density", "Genotype")

gen2.temp1 <- data.frame(gen2$radius_1, gen2$radius_2, gen2$radius_3, gen2$radius_4)
gen2.temp2 <- as.data.frame(unlist(gen2.temp1))
gen2.temp3 <- as.data.frame(unlist(data.frame(gen2$density_1, gen2$density_2, gen2$density_3, gen2$density_4)))
gen2.temp3$genotype <- c("588711")
gen2.fin <- cbind(gen2.temp2, gen2.temp3)
colnames(gen2.fin) <- c("Radius", "Density", "Genotype")

data <- rbind(gen1.fin, gen2.fin)
#data$Genotype <- gsub("588710", "Little", data$Genotype)
#data$Genotype <- gsub("588711", "Big", data$Genotype)
data$Genotype <- as.factor(data$Genotype)

# Plot density as box and whisker plot
gen1.den <- data[data$Genotype=="588710",]
gen1.counts <- as.data.frame(table(gen1.den$Density))
gen1.counts$Genotype <- c("588710")
gen2.den <- data[data$Genotype=="588711",]
gen2.counts <- as.data.frame(table(gen2.den$Density))
gen2.counts$Var1<-factor(gen2.counts$Var1,levels=c(0,1,3,5,7,9))
gen2.counts$Genotype <- c("588711")
gen2.counts<-rbind(gen2.counts,data.frame(Var1=c(0,1),Freq=0,Genotype=588711))
gen2.counts<-gen2.counts[order(gen2.counts$Var1),]

counts <- rbind(gen1.counts, gen2.counts)
colnames(counts) <- c("Density", "Counts", "Genotype")
#counts$Genotype <- gsub("588710", "Little", counts$Genotype)
#counts$Genotype <- gsub("588711", "Big", counts$Genotype)

########################## PLOT PHENOTYPE DATA ##########################

col.711 <- c("#024F4A")
col.710 <- c("#05B384")

#Plotting density as a a bar plot
# b <- ggplot(counts, aes(x=Density, y=Counts, fill=Genotype)) + 
#   geom_bar(stat='identity', position = 'dodge') +
#   theme_classic() +
#   scale_fill_manual(values = c(col.710, col.711), limits = c('588710', '588711')) +
#   theme(axis.title.x=element_text(colour="black", size=14),
#         axis.title.y=element_text(colour="black", size=14, margin = margin(t = 0, r = 10, b = 0, l = 0)),
#         axis.text.x=element_text(colour="black", size=14),
#         axis.text.y = element_text(colour="black", size=14),
#         legend.position = "right",
#         legend.text = element_text(size = 14),
#         legend.title = element_text(size=14)) +
#   scale_y_continuous(expand = c(0,0))

#Plotting density as a violiin plot
b <- ggplot(data, aes(x=factor(Genotype, levels = c('588710', '588711')), y=Density, fill=Genotype)) + 
  geom_violin(alpha=0.75) + geom_jitter(shape=16, position=position_jitter(0.2)) +
  theme_classic() +
  scale_fill_manual(values = c(col.710, col.711), limits = c('588710', '588711')) +
  theme(axis.title.x=element_text(colour="black", size=14),
        axis.title.y=element_text(colour="black", size=14, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.x=element_text(colour="black", size=14),
        axis.text.y = element_text(colour="black", size=14),
        legend.position = "none") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 10.5), breaks = c(0, 1, 3, 5, 7, 9) ) +
  ylab("Domatia Density") +
  xlab("Genotype") +
  geom_signif(
    comparisons = list(c('588710', '588711')),
    test= t.test,
    map_signif_level = TRUE, textsize = 6,
    y_position = 9.5) +
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 1,
               colour = "black")

# Plot radius as a violin plot
a <- ggplot(data, aes(x=factor(Genotype, levels = c('588710', '588711')), y=Radius, fill=Genotype)) + 
  geom_violin(alpha=0.75) + geom_jitter(shape=16, position=position_jitter(0.2)) +
  theme_classic() +
  scale_fill_manual(values = c(col.710, col.711), limits = c('588710', '588711')) +
  theme(axis.title.x=element_text(colour="black", size=14),
        axis.title.y=element_text(colour="black", size=14, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.x=element_text(colour="black", size=14),
        axis.text.y = element_text(colour="black", size=14),
        legend.position = "none") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 2.8), breaks = c(0, 0.5, 1.0, 1.5, 2.0, 2.5) ) +
  ylab("Domatia Radius (mm)") +
  xlab("Genotype") +
  geom_signif(
    comparisons = list(c('588710', '588711')),
    test= t.test,
    map_signif_level = TRUE, textsize = 6,
    y_position = 2.5) +
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 1,
               colour = "black")
  
plot_grid(a, NULL, b, NULL, nrow=1, ncol=4, rel_widths = c(1,0.125,1, 0.1))

########################## STATS ##########################
# Run t-test on radius data
a = data[data$Genotype=="588710",1]
b = data[data$Genotype=="588711",1]
ttest1 <- t.test(a,b)

# Run Mann-Whitney test on density data (which are discrete and not normally distributed)
c = data[data$Genotype=="588710",2]
d = data[data$Genotype=="588711",2]
wilcox1 <- wilcox.test(c, d, alternative = "less")
