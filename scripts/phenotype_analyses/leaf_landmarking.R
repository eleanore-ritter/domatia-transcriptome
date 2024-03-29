# Set working directory and load packages
#setwd("C:/Users/rittere5/OneDrive - Michigan State University/Vitis-domatia/leaf-landmarking/") #Work working directory
setwd("C:/Users/elean/OneDrive - Michigan State University/Vitis-domatia/leaf-landmarking/")
library(ggplot2)
library(shapes)
library(reshape2)
library(dplyr)
library(gridExtra)
library("cowplot")
library(plotrix)
library(stringr)
library(ggsignif)
library(Momocs)
library(ggnewscale)
library(viridis)

########################## 1 LOAD DATA IN AND RUN INITIAL ANALYSES ##########################

#Load data in
files<-list.files()[grepl('.txt',list.files())]
tmp1<-read.csv(files[1],sep="",header=FALSE)
coords<-array(dim=c(dim(tmp1),length(files)))
for(i in seq_along(files)){
  coords[,,i]<-as.matrix(read.csv(files[i],sep="",header=FALSE))
}

#Run procGPA and make eigenleaves with shapepca
gpa <- procGPA(coords, reflect=TRUE)
scaled <- procGPA(coords, reflect=TRUE, scale = TRUE)
shapepca(gpa, joinline=c(1:21,1)) #Need to fix joinline
shapepca(gpa, pcno=c(1,2,3), joinline=c(1,6, 14, 15, 16, 17, 18, 19, 20, 21, 13, 4), mag=0.5)

#Extract out and plot rotated points
tmp2 <- gpa$rotated

rotated <- data.frame(id=rep(files,each=dim(tmp2)[1]),
           landmark=seq_len(dim(tmp2)[1]),
           do.call(rbind,asplit(tmp2,3)))

ggplot(rotated, aes(x=X1, y=X2)) +
  geom_point(size=2, aes(color=landmark)) +
  theme_classic() +
  coord_fixed() +
  theme(aspect.ratio = 1)

########################## 2 CHECK ONE, SPECIFIC FILE IF NEEDED ##########################
#Load in file
c1 <- read.csv("588710-1-A.txt", sep="" ,header = FALSE)
#Visualize points
ggplot(c1, aes(x=V1, y=V2)) +
  geom_point(size=2) +
  theme_classic() +
  coord_fixed() +
  theme(aspect.ratio = 1)

########################## 3 TEST MEAN SHAPE DIFFERENCES BETWEEN SAMPLES ##########################
# Set working directory and load packages
#setwd("C:/Users/rittere5/OneDrive - Michigan State University/Vitis-domatia/leaf-landmarking/") #Work working directory
setwd("C:/Users/elean/OneDrive - Michigan State University/Vitis-domatia/leaf-landmarking/")

#Rerun procGPA for each sample
## s710
s710.files<-list.files()[grepl('588710',list.files())]
s710.tmp1<-read.csv(s710.files[1],sep="",header=FALSE)
s710.coords<-array(dim=c(dim(s710.tmp1),length(s710.files)))
for(i in seq_along(s710.files)){
  s710.coords[,,i]<-as.matrix(read.csv(s710.files[i],sep="",header=FALSE))
}
s710.gpa <- procGPA(s710.coords, reflect=TRUE, scale = TRUE)

## s711
s711.files<-list.files()[grepl('588711',list.files())]
s711.tmp1<-read.csv(s711.files[1],sep="",header=FALSE)
s711.coords<-array(dim=c(dim(s711.tmp1),length(s711.files)))
for(i in seq_along(s711.files)){
  s711.coords[,,i]<-as.matrix(read.csv(s711.files[i],sep="",header=FALSE))
}
s711.gpa <- procGPA(s711.coords, reflect=TRUE, scale = TRUE)

#Test mean shape difference between samples
tms <- testmeanshapes(s710.coords, s711.coords, scale = TRUE)

# Plot mean shapes together
#Run procGPA once to get points scaled and rotated the same
#Load data in
files<-list.files()[grepl('.txt',list.files())]
tmp1<-read.csv(files[1],sep="",header=FALSE)
coords<-array(dim=c(dim(tmp1),length(files)))
for(i in seq_along(files)){
  coords[,,i]<-as.matrix(read.csv(files[i],sep="",header=FALSE))
}

#Run procGPA and make eigenleaves with shapepca
gpa.scaled <- procGPA(coords, reflect=TRUE, scale = TRUE)
scaled <- gpa.scaled$rotated

#Extract scaled and rotated points by sample
which(grepl('588710', as.character(files)))
s710.scaled <- scaled[1:21,1:2,1:25] 
s710.gpa.scaled <- procGPA(s710.scaled, reflect=TRUE, scale = TRUE)

which(grepl('588711', as.character(files)))
s711.scaled <- scaled[1:21,1:2,26:54] 
s711.gpa.scaled <- procGPA(s711.scaled, reflect=TRUE, scale = TRUE)

# Get mean shapes
s710.mshape <- as.data.frame(s710.gpa.scaled$mshape)
s710.mshape$sample <- c("588710")
s711.mshape <- as.data.frame(s711.gpa.scaled$mshape)
s711.mshape$sample <- c("588711")

#Order data properly before combining
#Make vectors for geom_path() order
orderout <- c("1", "6", "14", "15", "16", "17", "18", "19", "20", "21", "13", "4") #Make vector for order of landmarks to join for leaf edge
orderin <- c("1", "6", "14", "5", "15", "7", "2", "9", "17", "8", "18",
             "10", "3", "12", "20", "11", "21", "13", "4") #Make vector for order of landmarks to join for veins

#Plot s710 and color by leaf node
s710.mshape$landmark <- rownames(s710.mshape)
s710.mshape$landmark <- as.character(s710.mshape$landmark)#Reorder nodes to be in proper joinline order

#Make dataframe with points arranged for outer line
s710.mshape.O <- s710.mshape %>% arrange(factor(landmark, levels = orderout))
s710.mshape.O <- s710.mshape.O[s710.mshape.O$landmark %in% orderout, ]
s710.mshape.O$line <- c("outer")

#Make dataframe with points arranged for vein lines
s710.mshape.I <- s710.mshape %>% arrange(factor(landmark, levels = orderin))
s710.mshape.I <- s710.mshape.I[s710.mshape.I$landmark %in% orderin, ]
s710.mshape.I$line <- c("inner")

#Merge two dataframes
s710.mshape.A <- rbind(s710.mshape.I, s710.mshape.O)

#Plot s711 and color by leaf node
s711.mshape$landmark <- rownames(s711.mshape)
s711.mshape$landmark <- as.character(s711.mshape$landmark)#Reorder nodes to be in proper joinline order

#Make dataframe with points arranged for outer line
s711.mshape.O <- s711.mshape %>% arrange(factor(landmark, levels = orderout))
s711.mshape.O <- s711.mshape.O[s711.mshape.O$landmark %in% orderout, ]
s711.mshape.O$line <- c("outer")

#Make dataframe with points arranged for vein lines
s711.mshape.I <- s711.mshape %>% arrange(factor(landmark, levels = orderin))
s711.mshape.I <- s711.mshape.I[s711.mshape.I$landmark %in% orderin, ]
s711.mshape.I$line <- c("inner")

#Merge two dataframes
s711.mshape.A <- rbind(s711.mshape.I, s711.mshape.O)


# Plot mean shapes
# Merge dataframes by variety
both.mshape <- rbind(s710.mshape.A, s711.mshape.A)

#Set colors
col1 <- c("#05B384")
col2 <- c("#024F4A")

# Plot mean shapes
p.both.mshape <- ggplot() +
  geom_path(s710.mshape.A, mapping = aes(x=V1, y=V2, group=line, color='588710'), alpha=0.75, size = 1.5) +
  geom_path(s711.mshape.A, mapping = aes(x=V1, y=V2, group=line, color='588711'), alpha=0.75, size = 1.5) +
  theme_classic() +
  coord_fixed() +
  scale_color_manual(values = c(col1, col2)) +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        plot.margin=unit(c(-1,-2,-1,-2), "pt"),
        legend.position = "bottom") +
  guides(col = guide_legend(override.aes = list(alpha = 1)))

########################## 4 PLOTTING PCA FOR SHAPE DATA ##########################
# Set working directory and load packages
#setwd("C:/Users/rittere5/OneDrive - Michigan State University/Vitis-domatia/leaf-landmarking/") #Work working directory
setwd("C:/Users/elean/OneDrive - Michigan State University/Vitis-domatia/leaf-landmarking/")

#Set colors
col1 <- c("#05B384")
col2 <- c("#024F4A")

#Load data in by genotype
m.files<-list.files()[grepl('58871',list.files())]
m.tmp1<-read.csv(m.files[1],sep="",header=FALSE)
m.coords<-array(dim=c(dim(m.tmp1),length(m.files)))
for(i in seq_along(m.files)){
  m.coords[,,i]<-as.matrix(read.csv(m.files[i],sep="",header=FALSE))
}

#Run procGPA and restructure data
m.gpa <- procGPA(m.coords, reflect=TRUE, scale = TRUE)
shapepca(m.gpa, pcno=c(1,2,3,4), joinline=c(1,6, 14, 15, 16, 17, 18, 19, 20, 21, 13, 4), mag=0.5, type = "r") #c equals 3 when type="r"
m.pc <- as.data.frame(m.gpa$scores)
m.sample <- gsub( "-.*", "", m.files)
m.pc$sample <- m.sample
mshoot <- substr(m.files, 8, 8)
mtmpnode1 <- gsub( ".*-", "", m.files)
myLetters <- LETTERS[1:26]
mtmpnode2 <- gsub(".txt", "", mtmpnode1)
mtmpnode3 <- match(mtmpnode2, myLetters)
mnode <- as.numeric(mtmpnode3)
rm(mtmpnode1, mtmpnode2, mtmpnode3)
m.pc$shoot <- mshoot
m.pc$node <- mnode

#Separate datasets for plotting
g710pc <-m.pc[m.pc$sample=="588710",]
g711pc <-m.pc[m.pc$sample=="588711",]

#Plot both genotypes samples
mpcA <- ggplot(g710pc, aes(x=PC1, y=PC2)) +
  geom_point(color=col1, size=2, alpha = 0.75) +
  stat_ellipse(data=g710pc,aes(x=PC1, y=PC2),type = "norm", color=col1, size=1) +
  geom_point(data=g711pc, aes(x=PC1, y=PC2), color=col2, size=2, alpha = 0.75, show.legend = F) +
  stat_ellipse(data=g711pc,aes(x=PC1, y=PC2),type = "norm", color=col2, size=1) +
  theme_classic() +
  labs(x = "PC1: 41.2%", y = "PC2: 18.3%") +
  theme(axis.title.x=element_text(colour="black", size=12),
        axis.title.y=element_text(colour="black", size=12),
        axis.text.x=element_text(colour="black", size=10),
        axis.text.y = element_text(colour="black", size=10),
        legend.position = "bottom",
        legend.justification = "center",
        legend.title = element_text(colour="black", size=10),
        legend.key.size = unit(0.25, "cm"),
        legend.text = element_text(colour="black", size=10))

mpcB <- ggplot(g710pc, aes(x=PC3, y=PC4)) +
  geom_point(color=col1, size=2, alpha = 0.75) +
  stat_ellipse(data=g710pc,aes(x=PC3, y=PC4),type = "norm", color=col1, size=1) +
  geom_point(data=g711pc, aes(x=PC3, y=PC4), color=col2, size=2,alpha = 0.75, show.legend = F) +
  stat_ellipse(data=g711pc,aes(x=PC3, y=PC4),type = "norm", color=col2, size=1) +
  theme_classic() +
  labs(x = "PC3: 9.3%", y = "PC4: 6.7%") +
  theme(axis.title.x=element_text(colour="black", size=12),
        axis.title.y=element_text(colour="black", size=12),
        axis.text.x=element_text(colour="black", size=10),
        axis.text.y = element_text(colour="black", size=10),
        legend.position = "bottom",
        legend.justification = "center",
        legend.title = element_text(colour="black", size=10),
        legend.key.size = unit(0.25, "cm"),
        legend.text = element_text(colour="black", size=10))

plot_grid(mpcA, NULL, mpcB, ncol=1, nrow=3, rel_heights = c(1,0.1, 1), align = "v", labels = "auto")

pc.grid <- plot_grid(mpcA, NULL, mpcB, ncol=1, nrow=3,  rel_heights = c(1,0.1, 1), align = "v", labels = c('b', '', 'c'))

plot_grid(p.both.mshape, pc.grid, ncol=2, nrow=1, align = "tblr", rel_widths = c(1,.85), labels = c('a', '', '')) 