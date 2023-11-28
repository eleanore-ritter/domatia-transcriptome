# Set working directory and load necessary packages
setwd("C:/Users/rittere5/OneDrive - Michigan State University/Vitis-domatia/")
#setwd("C:/Users/elean/OneDrive - Michigan State University/Vitis-domatia/")

library(dplyr)
library(tidyr)
library(ggplot2)

# Read in differentially expressed genes
dom710 <- read.csv("differentially_expressed_genes/DE_genes_Domatia_V_Leaf_588710.csv")
dom711 <- read.csv("differentially_expressed_genes/DE_genes_Domatia_V_Leaf_588711.csv")

# Read in Arabidopsis genes based on GO terms
cellwall710 <- read.csv("Nov-2023-GO-terms/588710-upregulated-cell-wall-biosyn-modif.txt",
                        "\t", skip = 1, header = FALSE)
all.cellwall710 <- merge(cellwall710, dom710, by.x="V2",by.y="ensembl_gene_id")

cellwall711 <- read.csv("Nov-2023-GO-terms/588711-upregulated-cell-wall-biosyn-modif.txt",
                        "\t", skip = 1, header = FALSE)
all.cellwall711 <- merge(cellwall711, dom711, by.x="V2",by.y="ensembl_gene_id")

xylan710 <- read.csv("Nov-2023-GO-terms/588710-upregulated-xylan.txt",
                     "\t", skip = 1, header = FALSE)
all.xylan710 <- merge(xylan710, dom710, by.x="V2",by.y="ensembl_gene_id")

xylan711 <- read.csv("Nov-2023-GO-terms/588711-upregulated-xylan.txt",
                     "\t", skip = 1, header = FALSE)
all.xylan711 <- merge(xylan711, dom711, by.x="V2",by.y="ensembl_gene_id")

pectin710 <- read.csv("Nov-2023-GO-terms/588710-upregulated-pectin.txt",
                     "\t", skip = 1, header = FALSE)
all.pectin710 <- merge(pectin710, dom710, by.x="V2",by.y="ensembl_gene_id")

pectin711 <- read.csv("Nov-2023-GO-terms/588711-upregulated-pectin.txt",
                      "\t", skip = 1, header = FALSE)
all.pectin711 <- merge(pectin711, dom711, by.x="V2",by.y="ensembl_gene_id")

lignin710 <- read.csv("Nov-2023-GO-terms/588710-upregulated-lignin.txt",
                      "\t", skip = 1, header = FALSE)
all.lignin710 <- merge(lignin710, dom710, by.x="V2",by.y="ensembl_gene_id")

lignin711 <- read.csv("Nov-2023-GO-terms/588711-upregulated-lignin.txt",
                      "\t", skip = 1, header = FALSE)
all.lignin711 <- merge(lignin711, dom711, by.x="V2",by.y="ensembl_gene_id")

# Read in Arabidopsis metabolism genes based on GO terms
aa710 <- read.csv("GO-term-enrichment/for-heatmaps/metabolism/588710-amino-acid.txt",
                  "\t", skip = 1, header = FALSE)
all.aa710 <- merge(aa710, dom710, by.x="V2",by.y="ensembl_gene_id")

aa711 <- read.csv("GO-term-enrichment/for-heatmaps/metabolism/588711-amino-acid.txt",
                  "\t", skip = 1, header = FALSE)
all.aa711 <- merge(aa711, dom711, by.x="V2",by.y="ensembl_gene_id")
