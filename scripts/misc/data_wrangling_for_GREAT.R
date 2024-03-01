# Set working directory and load necessary packages
#setwd("C:/Users/rittere5/OneDrive - Michigan State University/Vitis-domatia/")
setwd("C:/Users/elean/OneDrive - Michigan State University/Vitis-domatia/")

library("plotly")

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
flower <- c(
#Seems highly expressed in stem   "Vitvi02g00060",
            "Vitvi01g01028",
            "Vitvi04g01440",
            "Vitvi06g00463",
            "Vitvi04g02117",
            "Vitvi08g01030",
#Seems highly expressed in stem            "Vitvi13g00740",
            "Vitvi15g00520",
            "Vitvi02g00512",
            "Vitvi04g02123",
            "Vitvi07g01906",
            "Vitvi18g03417",
            "Vitvi19g00264",
# Interesting gene, but seems to be moderately expressed in leaves            "Vitvi11g00247",
            "Vitvi11g00231",
            "Vitvi10g00854",
            "Vitvi10g00644",
            "Vitvi13g00829",
            "Vitvi12g02361",
            "Vitvi12g00725",
            "Vitvi15g00776"
# Highly expressed in Arabidopsis flowers, but somwhat expressed in Vvi stem and leaves            ,"Vitvi14g01928"
            )

flower.vvi <- as.matrix(counts[row.names(counts) %in% flower,])

temp1 <- as.data.frame(merge(dict, dom, by.x = "V1", by.y = "Gene")) #Good for searching gene IDs
flower <- gsub("$", ".t01", flower)
flower <- as.data.frame(flower)
flower.genes <- temp1[temp1$V2 %in% flower$flower,]

heatmap(flower.vvi, Rowv = NA, Colv = NA, margins=c(10,10))

p <- plot_ly(x=colnames(flower.vvi), y=rownames(flower.vvi), z = flower.vvi, type = "heatmap")
p

# Log transform
flower.vvi.log2 <- log2(flower.vvi + 1)
flower.vvi.Z <- t(scale(t(flower.vvi.log2)))
col.order <- c(
#Floral  
  "Vitvi15g00776",
  "Vitvi10g00854",
  "Vitvi08g01030",
  "Vitvi01g01028",
  "Vitvi11g00231",
  "Vitvi15g00520",
# JMT, BSMT, or AAP2 (defense-related)
  "Vitvi07g01906",
  "Vitvi04g02123",
  "Vitvi04g02117",
  "Vitvi12g00725",
#Trichome  
  "Vitvi04g01440",
#Metabolism  
  "Vitvi02g00512",
  "Vitvi19g00264",
#Misc
  "Vitvi13g00829",
  "Vitvi06g00463",
  "Vitvi10g00644",
  "Vitvi12g02361"
  
)
flower.vvi.Z <- flower.vvi.Z[col.order ,]
heatmap(flower.vvi.Z, Rowv = NA, Colv = NA, margins=c(10,10))

.lin.interp<-function(x,length.out){
  xx<-seq(1,length(x),length.out=length.out)
  in.inds<-xx%in%(1:length(xx))
  out<-rep(NA,length.out)
  out[in.inds]<-x[xx[in.inds]]
  xx<-xx[!in.inds]
  out[!in.inds]<-(x[ceiling(xx)]-x[floor(xx)])/(ceiling(xx)-floor(xx))*(xx-floor(xx))+x[floor(xx)]
  out
}

alter.cols<-function(x,alpha=NA,mix=NA,wgt=0.5){
  lens<-lengths(list(x,alpha,mix,wgt))
  #prepare mix
  if(any(!is.na(mix))){
    mix.len<-max(lens[-c(1,2)])
    out.len<-max(lens)
    tmp<-rep(mix,length.out=mix.len)
    wgt[wgt<0]<-0
    wgt[wgt>1]<-1
    wgt<-rep(wgt,length.out=mix.len)
    wgt[is.na(tmp)]<-0
    mix<-col2rgb(mix,alpha=TRUE)[,rep(seq_len(lens[3]),length.out=mix.len),drop=FALSE]
    mix<-sweep(mix,2,wgt,'*')
    mix<-mix[,rep(seq_len(mix.len),length.out=out.len),drop=FALSE]
  }else{
    out.len<-max(lens[c(1,2)])
    mix<-NULL
  }
  #prepare x
  x<-col2rgb(x,alpha=TRUE)[,rep(seq_len(lens[1]),length.out=out.len),drop=FALSE]
  if(!is.null(mix)){
    x<-sweep(x,2,1-wgt,'*')
    x<-x+mix
  }
  #prepare alpha
  if(any(!is.na(alpha))){
    pre<-x[4,,drop=FALSE]
    alpha[alpha<0]<-0
    alpha[alpha>1]<-1
    x[4,]<-255*alpha
    nas<-is.na(x[4,,drop=FALSE])
    x[4,nas]<-pre[nas]
  }
  do.call(rgb,c(asplit(x,1),maxColorValue=255))
}

legend.evorates<-function(mat,location=c('bottomleft','topleft','bottomright','topright'),
                          exp=FALSE,exp.txt=FALSE,
                          col=hcl.colors(12,"YlOrRd",rev=TRUE),val.range=NULL,res=100,
                          alpha=NA,breaks=NULL,select.levels=NULL,
                          box.dims=NULL,box.offset=NULL,box.scale=1,
                          txt.col=NULL,main=NULL,...){
  if(exp){
    exp.txt<-FALSE
  }
  try.location<-try(match.arg(location,c('bottomleft','topleft','bottomright','topright')),silent=T)
  if(inherits(try.location,'try-error')){
    stop(location," is not an available named position to put the legend: please specify one of the following: 'bottomleft', 'topleft', 'bottomright', or 'topright'")
  }
  location<-try.location
  gen.args<-graphics:::.Pars
  poly.args<-names(formals(polygon))
  poly.args<-poly.args[-which(poly.args=='...')]
  text.args<-names(formals(text.default))
  text.args<-text.args[-which(text.args=='...')]
  if(is.null(breaks)){
    if(is.null(val.range)){
      val.range<-range(mat,na.rm=TRUE)
    }
    colramp<-colorRampPalette(col,alpha=T)(res)
    colramp<-alter.cols(colramp,alpha=.lin.interp(alpha,length(colramp)))
  }else{
    colramp<-colorRampPalette(col,alpha=T)(length(breaks)+1)
    colramp<-alter.cols(colramp,alpha=.lin.interp(alpha,length(colramp)))
    if(is.null(select.levels)){
      select.levels<-1:length(colramp)
    }
    select.levels<-select.levels[select.levels>=1&select.levels<=(length(breaks)+1)]
    select.levels<-sort(select.levels)
    colramp<-colramp[select.levels]
  }
  bds<-par('usr')
  bds.dims<-c(diff(bds[1:2]),diff(bds[3:4]))
  if(length(box.dims)==0){
    box.dims<-rep(NA,2)
  }else if(length(box.dims)==1){
    box.dims<-c(box.dims,NA)
  }
  box.dims<-box.scale*ifelse(is.na(box.dims),bds.dims/c(30,5),box.dims)
  if(length(box.offset)==0){
    box.offset<-rep(NA,2)
  }else if(length(box.offset)==1){
    box.offset<-c(box.offset,NA)
  }
  box.offset<-ifelse(is.na(box.offset),bds.dims/c(8,20),box.offset)
  x.offset<-box.offset[1]
  y.offset<-box.offset[2]
  if(grepl('right',location)){
    x.offset<- bds.dims[1]-box.dims[1]-box.offset[1]
  }
  if(grepl('top',location)){
    y.offset<- bds.dims[2]-box.dims[2]-box.offset[2]
  }
  coords<-list(x=c(0,box.dims[1],box.dims[1],0)+x.offset+bds[1],y=c(0,0,box.dims[2],box.dims[2])+y.offset+bds[3])
  # box.midpt<-sapply(coords,mean)
  # coords$x<-(coords$x-box.midpt[1])*box.scale+box.midpt[1]
  # coords$y<-(coords$y-box.midpt[2])*box.scale+box.midpt[2]
  y.int<-seq(coords$y[2],coords$y[3],length.out=length(colramp)+1)
  for(i in 1:length(colramp)){
    do.call(polygon,
            c(x=list(coords$x),
              y=list(c(y.int[i],y.int[i],y.int[i+1],y.int[i+1])),
              border=NA,
              col=colramp[i],
              list(...)[!(names(list(...))%in%c('x','y','border','col','adj'))&
                          names(list(...))%in%gen.args]))
  }
  do.call(polygon,
          c(x=list(coords$x),
            y=list(coords$y),
            col=NA,
            list(...)[!(names(list(...))%in%c('x','y','col','adj'))&
                        names(list(...))%in%c(gen.args,poly.args)]))
  side<-NA
  if(hasArg(side)){
    if(list(...)$side<=2&list(...)$side>=1){
      side<-list(...)$side
    }
  }
  if(is.na(side)){
    side<-if(grepl('left',location)) 2 else 1
  }
  txt.args<-list(...)[!(names(list(...))%in%c('x','y','labels'))&names(list(...))%in%c(gen.args,text.args)]
  if(is.null(txt.args$adj)&is.null(txt.args$pos)){
    txt.args$pos<-side*2
  }
  if(!is.null(txt.col)){
    txt.args$col<-txt.col
  }
  if(is.null(breaks)){
    if(val.range[2]-val.range[1]==0){
      labels<-val.range[1]
    }else{
      labels<-pretty(seq(val.range[1],val.range[2],length.out=100))
      labels<-labels[2:(length(labels)-1)]
    }
    y.pos<-coords$y[2]+(labels-val.range[1])/
      (diff(val.range))*
      (coords$y[3]-coords$y[2])
    if(exp.txt){
      labels<-format(exp(labels),digits=1)
    }
    do.call(text,
            c(x=coords$x[side],
              y=list(y.pos),
              labels=list(labels),
              txt.args))
  }else{
    labels<-paste(signif(breaks[-length(breaks)],3),signif(breaks[-1],3),sep=' - ')
    labels<-c(paste('<',signif(breaks[1],3)),labels,paste('>',signif(breaks[length(breaks)],3)))
    labels<-labels[select.levels]
    y.pos<-apply(cbind(y.int[-length(y.int)],y.int[-1]),1,mean)
    do.call(text,
            c(x=coords$x[side],
              y=list(y.pos),
              labels=list(labels),
              txt.args))
  }
  if(is.null(main)){
    main=""
  }
  main.side<-NA
  if(hasArg(main.side)){
    if(list(...)$main.side<=4&list(...)$main.side>=1){
      main.side<-list(...)$main.side
    }
  }
  if(is.na(main.side)){
    main.side<-if(grepl('left',location)) 2 else 4
  }
  main.args<-list(...)[!(names(list(...))%in%paste0('main.',c('x','y','labels','side')))&
                         names(list(...))%in%paste0('main.',c(gen.args,text.args))]
  names(main.args)<-gsub('main.','',names(main.args))
  txt.args<-txt.args[!names(txt.args)%in%c('srt','pos','adj','offset')]
  txt.args[names(main.args)%in%names(txt.args)]<-main.args[names(main.args)%in%names(txt.args)]
  main.args<-c(txt.args,main.args[!(names(main.args)%in%names(txt.args))])
  if(is.null(main.args$adj)&is.null(main.args$pos)){
    main.args$pos<-main.side
    if(main.side%in%c(2,4)&is.null(main.args$srt)){
      main.args$srt<-90
      main.args$pos<-NULL
      main.args$adj<-c(0.5,c(-0.3,1.3)[main.side/2])
    }
  }
  if(main.side%in%c(2,4)&is.null(main.args$srt)){
    main.args$srt<-90
  }
  if(is.null(main.args$cex)){
    main.args$cex<-1
  }
  if(main.side%in%c(2,4)){
    x.pos<-coords$x[main.side/2]
    y.pos<-mean(coords$y[2:3])
    if(main.side==4){
      x.pos<-x.pos+bds.dims[1]/50
    }else{
      x.pos<-x.pos-bds.dims[1]/50
    }
  }else{
    x.pos<-mean(coords$x[1:2])
    y.pos<-coords$y[main.side/2+1.5]
    if(main.side==3){
      y.pos<-y.pos+bds.dims[2]/50
    }else{
      y.pos<-y.pos-bds.dims[2]/50
    }
  }
  do.call(text,
          as.list(c(x=x.pos,
                    y=y.pos,
                    labels=list(main),
                    main.args)))
  invisible(coords)
}

legend.evorates(flower.vvi.Z,location="bottomright",box.offset=c(0,0))

#plot_ly(x=colnames(flower.vvi.Z), y = rownames(flower.vvi.Z), z = flower.vvi.Z, type = "heatmap")      

temp <- log2(genes.mat + 1)
genes.mat1 <- temp[300:350,]
plot_ly(x=colnames(genes.mat1), y = rownames(genes.mat1), z = genes.mat1, type = "heatmap") 

######################## HOW MANY DOMATIA GENES ARE EXPRESSED IN FLOWERS? ########################

countsv <- unlist(counts)
quantile(countsv,prob=0.25) # the 25% quantile of expressed reads is 17.45 reads

data1FL1 <- data1 %>% 
  filter(rowAny(across(where(is.numeric), ~ .x > 17.45)))
data1FL2 <- data1 %>% filter_all( ~ .x > 17.45)
