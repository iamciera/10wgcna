## WGCNA with bootstrap
## 4.06.2015.a
## Ciera Martinez

#I started with a copy of WGCNA.bootstrap.Ciera.03.27.15.c.  I am
#using Yasu's data. I want to make network with JUST cluster 35 and 17.
#I am testing if these two different areas cluster differently.
#Since each cluster should network differently. 
# [ ] Add visualization of the clusters at the end. 

#Cluster 17 is a mind fuck.  Something is wrong with it. 


## library
library(WGCNA)
options(stringsAsFactors  =  FALSE)
#enableWGCNAThreads()
ALLOW_WGCNAT_THREADS = 4
library (igraph)
library(ggplot2)
library(reshape)

library(rgl)
library(tcltk2)


## Dan Koenig 
counts <- read.delim("../data/GSE45774_rpkm_all.txt", header = TRUE)
colnames(counts)[1] <- "Gene_ID"
head(counts)
dim(counts)


#Read in lists of genes from both SOM and superSOM
#SOM_analysis9.5.csv is the large WT only analysis.

SOM <- read.csv("../../08SOM/lcmSOM/data/SOM_analysis9.5forNetwork.csv", header = TRUE) 
colnames(SOM)
head(SOM)
#Subset only the columns that specify clusters and genes
SOMsub <- SOM[,c(2,21,22)] 
names(SOMsub)


#Subsetted for the interesting genes

head(SOMsub)

cluster35 <- subset(SOMsub, som.unit.classif == "35")
dim(cluster35)


#isolate only gene names
SOMclusters <- as.data.frame(cluster35[,1])
colnames(SOMclusters)[1] <- "Gene_ID"
head(SOMclusters)

#remove duplicates
dim(SOMclusters)
SOMclusters <- unique(SOMclusters)
dim(SOMclusters) 
#no duplicates

#Now merge with Yasu's table to get only genes I am interested in.
#First rename 1st column in counts to gene for merging

#Merge
dim(counts)
head(counts)
dim(SOMclusters)

merged <- merge(SOMclusters, counts, by = "Gene_ID")

dim(merged)

#write.csv(countsUniq, file = "countsUniq.csv", row.names = FALSE)
#countsUniq <- read.csv("countsUniq.csv", row.name = 1) 

counts <- merged
counts[is.na(counts)] <- 0
str(counts)

dim(counts) 
head(counts)

genes <- counts[,1]

# transpose all but the first column (name)
counts.t <- as.data.frame(t(counts[,-1]))

colnames(counts.t) <- genes

head(counts.t)
dim(counts.t)

#counts.lt=t(log(counts+1)) 

############################################################
## Bootstrapping for hub gene prediction
B=100  ## select number of bootstrap resamples

powers  =  c(c(3:50)) #if you get error in the bootsrapping, you might need to the maximum value here.
result=matrix(nrow=ncol(counts.t), ncol=B)

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
qplot(average, sd, data=result.g)  #what should I be seeing here
colnames(result.g) <- c("ave.rank","sd.rank")

result.o <- result.g[order(result.g$ave.rank),]
top.hub <- rownames(result.o[1:400,]) # top  hub genes

# save
#save.image(file=paste("boot",B,".WGCNA.Rdata",sep=""))

## visualization
# Choose a set of soft-thresholding powers
powers  =  c(1:30)
#  Call  the  network  topology  analysis  function
sft=pickSoftThreshold(counts.t,powerVector=powers,RsquaredCut=0.9,verbose=5) #power must be between 1 and 30.

# create TOM (topological overlap matrix)
TOM =TOMsimilarityFromExpr(counts.t,power=14) # power=14 shows R^2=0.9
colnames(TOM)=genes
rownames(TOM)=genes

dim(TOM)
head(TOM,1)

# extract top hub genes, 
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
pdf("top.hub.test.04062015.pdf", useDingbats=FALSE) 
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

