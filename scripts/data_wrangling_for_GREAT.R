# Set working directory and load necessary packages
setwd("C:/Users/rittere5/OneDrive - Michigan State University/Vitis-domatia/")
#setwd("C:/Users/elean/OneDrive - Michigan State University/Vitis-domatia/")

######################## DATA WRANGLING FOR GREAT ########################

# Load data
dict <- read.csv("Vriparia-Vvinifera.tsv", sep='\t', header = FALSE)
df1 <- read.csv("differentially_expressed_genes/DE_genes_Domatia_V_Leaf_588710.csv")
row.names(df1) <- df1$Gene
df2 <- read.csv("differentially_expressed_genes/DE_genes_Domatia_V_Leaf_588711.csv")
row.names(df2) <- df2$Gene


# Shared between domatia and domatia
dom <- df1[rownames(df1)%in%rownames(df2),]

# Get Vvinifera orthologs to differentially expressed genes
vvi <- as.data.frame(merge(dict, dom, by.x = "V1", by.y = "Gene"))
vvi <- vvi$V2
vvi <- as.data.frame(gsub(".t01", "", vvi))
colnames(vvi) <- c("V1")
vvi<-t(vvi)

# Copy data for GREAT
writeClipboard(paste0(as.vector(vvi),collapse=","))

######################## DE ANALYSIS WITH GREAT DATA ########################
# Load GREAT data into R and reformat
data1 <- read.csv(file = "GREAT_analysis/GREAT-data/GREAT_RPKM_samples2023-07-27_inflorescence_flower.csv", row.names = "X")
data2 <- read.csv(file = "GREAT_analysis/GREAT-data/GREAT_RPKM_samples2023-07-27_leaf.csv", row.names = "X")
data3 <- read.csv(file = "GREAT_analysis/GREAT-data/GREAT_RPKM_samples2023-07-27_stem.csv", row.names = "X")

colnames(data1) <- paste("Inflorescence" , colnames(data1), sep=":")
colnames(data2) <- paste("Leaf" , colnames(data2), sep=":")
colnames(data3) <- paste("Stem" , colnames(data3), sep=":")

counts <- cbind(data1, data2, data3)
#counts$Gene <- row.names(counts)

genes.mat <- as.matrix(counts)
heatmap(genes.mat, Rowv = NA, Colv = NA,margins=c(10,10))

# List of genes that seem to be upregulated in inflorescences compared to leaf and stem tissue
flower <- c("Vitvi02g00512",
            "Vitvi04g00641",
            "Vitvi03g01477",
            "Vitvi04g02123",
            "Vitvi07g01906",
            "Vitvi07g00470",
            "Vitvi09g00463",
            "Vitvi18g03417",
            "Vitvi19g00369",
            "Vitvi19g00264",
            "Vitvi11g00247",
            "Vitvi11g00231",
            "Vitvi10g00854",
            "Vitvi10g00663",
            "Vitvi10g00644",
            "Vitvi13g00829",
            "Vitvi12g02361",
            "Vitvi12g00725",
            "Vitvi15g00776",
            "Vitvi14g01928",
            "Vitvi13g02416",
            "Vitvi18g02654")

flower.vvi <- as.matrix(counts[row.names(counts) %in% flower,])

temp1 <- as.data.frame(merge(dict, dom, by.x = "V1", by.y = "Gene")) #Good for searching gene IDs
flower <- gsub("$", ".t01", flower)
flower <- as.data.frame(flower)
flower.genes <- temp1[temp1$V2 %in% flower$flower,]

heatmap(flower.vvi, Rowv = NA, Colv = NA, margins=c(10,10))

p <- plot_ly(x=colnames(flower.vvi), y=rownames(flower.vvi), z = flower.vvi, type = "heatmap")
p
