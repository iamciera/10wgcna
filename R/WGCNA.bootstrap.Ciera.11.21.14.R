## WGCNA with bootstrap
## Nov 2014
## Ciera Martinez

## library
library(WGCNA)
options(stringsAsFactors  =  FALSE)
#enableWGCNAThreads()
ALLOW_WGCNAT_THREADS = 4
library (igraph)
library(ggplot2)
library(reshape)

#library(rgl)
#library(tcltk2)

## because WGCNA with soft threshold need > 12 samples, 

counts <- read.csv("../data/pnas.1402835111.sd01.csv", header = TRUE)
colnames(counts)


#Read in lists of genes from both SOM and superSOM

#######!!!!! You need to fix the gene names in the original.  They are truncated!
SOM <- read.csv("../../08SOM/lcmSOM/data/forNetwork/SOM_analysis9.5.csv", header = TRUE)
head(SOM)
superSOM <- read.csv("../../08SOM/lcmSOM/data/forNetwork/superSOM_analysis8.csv", header=TRUE)
head(superSOM)

#bring together lists from 

SOMgenes <- as.data.frame(SOM[,2])
colnames(SOMgenes)[1] <- "Gene_ID"
dim(SOMgenes)
duplicated(SOMgenes)

superSOMgenes <- as.data.frame(superSOM[,3])
dim(superSOMgenes)
colnames(superSOMgenes)[1] <- "Gene_ID"
duplicated(superSOMgenes)

SOMclusters <- rbind(superSOMgenes, SOMgenes)
dim(SOMclusters)
duplicated(SOMclusters)

#Merge
head(counts)
merged <- merge(SOMclusters, counts, by = "Gene_ID", all.x = TRUE)
dim(merged)
head(merged)

#remove rows that the have duplicate gene names.
 
dim(merged)
countsUniq <- unique(merged)
dim(countsUniq)

#Break point. Here are all the genes that are important in the gene clusters.
write.csv(countsUniq, file = "countsUniq.csv", row.names = FALSE)
countsUniq <- read.csv("countsUniq.csv", row.name = 1)

counts <- countsUniq 
counts[is.na(counts)] <- 0

#for later
genes <- rownames(counts) #setting the genes names, 

#transform data frame
counts.t <- t(counts)
str(counts.t)
head(counts.t)

#in Yasu's script, why is this used?
#counts[, c(2:64)] <- sapply(counts[, c(2:64)], as.numeric)
#counts.lt=t(log(counts+1)) 

############################################################
## Bootstrapping for hub gene prediction
B = 100  ## select number of bootstrap resamples

powers  =  c(c(3:50)) #if you get error in the bootsrapping, you might need to the maximum value here.
result = matrix(nrow=ncol(counts.t), ncol=B)

for (i in 1:B){
  
  set.seed(i*100+1)
  print(i)
  
  ##bootstrap resample
  sft.power=30
  
  while(sft.power>29 || is.na(sft.power)){#because TOM need power < 30 #softconnecity power < 14
    index.b=sample(x=1:nrow(counts.t), size=nrow(counts.t), replace=TRUE)
    Y.b=counts.t[index.b,]
    
    ##soft thresholding
    sft.b = pickSoftThreshold(Y.b,  powerVector=powers, RsquaredCut=0.9, verbose  =  5)
    sft.power = sft.b$powerEstimate
    
  }
  
  print(sft.power)
  ##TOM
  TOM.b = TOMsimilarityFromExpr(Y.b,power=sft.b$powerEstimate) #omega TOM-based connectivity
  hub.b = rowSums(TOM.b)
  #adj.b = adjacency(Y.b,power=sft.b$powerEstimate)
  #hub.b = rowSums(adj.b) #k connectivity
  
  result[,i]<-rank(-hub.b)
}
#There are NAs that are produced. Talk to Yasu about this.

#Annotation Files
annotation1<- read.delim("../../06diffGeneExp/analysis/data/ITAG2.3_all_Arabidopsis_ITAG_annotations.tsv", header=FALSE)  #Changed to the SGN human readable annotation
colnames(annotation1) <- c("ITAG", "SGN_annotation")
annotation2<- read.delim("../../06diffGeneExp/analysis/data/ITAG2.3_all_Arabidopsis_annotated.tsv")
annotation <- merge(annotation1,annotation2, by = "ITAG")
sub.annotation <- annotation[,c(1,4)] #Choose what you want to name them by

# #Merge with genes then call which column of genes you want.
table.genes <- as.data.frame(genes)
colnames(table.genes) <- "ITAG"

annotation.merge <- merge(table.genes,sub.annotation, by = "ITAG", all.x = TRUE)

#If NA, replace with ITAG column
annotation.merge$gene.name <- ifelse(is.na(annotation.merge$symbol), 
                                     annotation.merge$ITAG, 
                                     annotation.merge$symbol)

genes <-  annotation.merge$gene.name

#Make Ranking
row.names(result) <- genes
average <- rowMeans(result)
sd <- apply(result,1,function(d)sd(d))
result.n <- cbind(result,average,sd)
result.n <- as.data.frame(result.n)
result.g <- subset(result.n[,11:12])
qplot(average, sd, data=result.g) 
colnames(result.g) <- c("ave.rank","sd.rank")


result.o <- result.g[order(result.g$ave.rank),]
top.hub <- rownames(result.o[1:200,]) # top  hub genes

# save
save.image(file=paste("boot",B,".WGCNA.Rdata",sep=""))

## visualization
# Choose a set of soft-thresholding powers
powers  =  c(1:30)

#  Call  the  network  topology  analysis  function
sft=pickSoftThreshold(counts.t,powerVector=powers,RsquaredCut=0.9,verbose=5) #power must be between 1 and 30.

# create TOM (topological overlap matrix)
TOM =TOMsimilarityFromExpr(counts.t,power=sft$powerEstimate) # power=14 shows R^2=0.9
colnames(TOM)=genes
rownames(TOM)=genes

dim(TOM)
head(TOM)

# extract top hub genes, this is where the NAs are occuring. 
index.sub=is.element(genes, top.hub)
subTOM=TOM[index.sub,index.sub]

# only strong interaction is shown
h.subTOM = (subTOM>0.1)*subTOM # only > 0.1 TOM will be shown in network

subnet = graph.adjacency(h.subTOM,mode="undirected",weighted=TRUE,diag=FALSE)

between <- betweenness(subnet, normalized=TRUE)
#between.a<-merge(between,annotation, by.x="row.names", by.y="Sequence_name", all.x=T,sort=F)
#write.csv(between.a,"betweenness.csv")
head(between[order(-between)])

# visualization
V(subnet)$color <- "mediumturquoise"
#V(net)[community.fastgreedy$membership==1]$color <- "mediumturquoise"
v.label=rep("",length(V(subnet)))
v.label=V(subnet)$name
v.size=rep(3,length(V(subnet)))
V(subnet)$shape <- "circle"
pdf("top.hub.test.pdf", useDingbats=FALSE) 
plot(subnet, 
     layout=layout.graphopt, 
     vertex.size=v.size, 
     vertex.frame.color=NA,
     vertex.label=v.label, 
     vertex.label.cex=0.05,
     edge.color="gray57", 
     edge.width=E(subnet)$weight*0.1)
dev.off()


######################################

# Plot soft-thresholding powers
tiff("SoftThresholding.tif", width=16, height=8, unit="in",compression="lzw",res=100)
par(mfrow  =  c(1,2))
cex1  =  0.9
#  Scale-free  topology  fit  index  as  a  function  of  the  soft-thresholding  power
plot(sft$fitIndices[,1],  -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft  Threshold  (power)",ylab="Scale  Free  Topology  Model  Fit,signed  R^2",type="n",
     main  =  paste("Scale  independence"));
text(sft$fitIndices[,1],  -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
#  this  line  corresponds  to  using  an  R^2  cut-off  of  h
abline(h=0.90,col="red")
#  Mean  connectivity  as  a  function  of  the  soft-thresholding  power
plot(sft$fitIndices[,1],  sft$fitIndices[,5],
     xlab="Soft  Threshold  (power)",ylab="Mean  Connectivity",  type="n",
     main  =  paste("Mean  connectivity"))
text(sft$fitIndices[,1],  sft$fitIndices[,5],  labels=powers,  cex=cex1,col="red")
dev.off()

