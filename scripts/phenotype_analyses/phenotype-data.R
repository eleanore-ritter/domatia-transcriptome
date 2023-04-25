# Set working directory and load necessary packages
#setwd("C:/Users/rittere5/OneDrive - Michigan State University/Vitis-domatia/")
setwd("C:/Users/elean/OneDrive - Michigan State University/Vitis-domatia/")

library(ggplot2)
library(cowplot)
library(ggsignif)
library(tidyr)

# Load raw data
raw.data <- read.csv("domatia_trichomes_2022_final.csv")

# Extract out our two genotypes
gen1 <- raw.data[grep("588710",raw.data$plant_code),]
gen2 <- raw.data[grep("588711",raw.data$plant_code),]
data <- rbind(gen1, gen2)

# Add column for genotype
data$genotype <- data$plant_code
data$genotype <- gsub("-.*", "", data$genotype)

# Plot density_1 and radius_1 against dry mass to see if anything stands out
col.710 <- c("#E63946")
col.711 <- c("#457B9D")

p1 <- ggplot(data, aes(x = density_1, y = dry_mass..mg.)) +
  geom_point(aes(color = genotype), size=2, shape=19) +
  scale_color_manual(values = c(col.710, col.711)) +
  theme_classic() +
  labs(x = "Density 1", y = "Dry Mass (mg)") +
  theme(axis.title.x=element_text(colour="black", size=12),
        axis.title.y=element_text(colour="black", size=12),
        axis.text.x=element_text(colour="black", size=10),
        axis.text.y = element_text(colour="black", size=10),
        legend.position = "none")

p2 <- ggplot(data, aes(x = radius_1, y = dry_mass..mg.)) +
  geom_point(aes(color = genotype), size = 2, shape=19) +
  scale_color_manual(values = c(col.710, col.711)) +
  theme_classic()+
  labs(x = "Radius 1", y = "Dry Mass (mg)") +
  theme(axis.title.x=element_text(colour="black", size=12),
        axis.title.y=element_text(colour="black", size=12),
        axis.text.x=element_text(colour="black", size=10),
        axis.text.y = element_text(colour="black", size=10),
        legend.position = "right",
        legend.justification = "center",
        legend.title = element_blank(),
        legend.key.size = unit(0.25, "cm"),
        legend.text = element_text(colour="black", size=12))

plot_grid(p1, p2, nrow=1, ncol=2, rel_widths = c(1.6,2))

# ## Run t-test
# a = na.omit(data[data$genotype==588710,])
# b = na.omit(data[data$genotype==588711,])
# 
# ttest <- t.test(a,b) #Need to fix this, a ttest is probably not best!


