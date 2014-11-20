## WGCNA with bootstrap
## 2014 1 21

## library
library(WGCNA)
options(stringsAsFactors  =  FALSE)
enableWGCNAThreads()
library (igraph)
library(ggplot2)
library(reshape)

## Read in data: dodder cluster 5 with indicidual biological replicates
## because WGCNA with soft threshold need > 12 samples

counts <- read.csv("counts.csv",row.names=1)
dim(counts) #1498   40
head(counts)

genes=rownames(counts)
counts.t <- t(counts)
#counts.lt=t(log(counts+1))


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

row.names(result) <- genes
average <- rowMeans(result)
sd <- apply(result,1,function(d)sd(d))
result.n <- cbind(result,average,sd)
result.n <- as.data.frame(result.n)
result.g <- subset(result.n[,11:12])
qplot(average, sd, data=result.g) 
colnames(result.g) <- c("ave.rank","sd.rank")

# gene annotation
#annotation <- read.csv("Dodder_cdhitest95_TAIR_Blast2GO_annotation.csv", header=TRUE)
#result.a<-merge(result.g,annotation, by.x="row.names", by.y="Sequence_name", all.x=T,sort=F)
#write.csv(result.a,paste("boot",B,".WGCNA.hub.csv",sep=""))

result.o <- result.g[order(result.g$ave.rank),]
top.hub <- rownames(result.o[1:100,]) # top 100 hub genes

# save
save.image(file=paste("boot",B,".WGCNA.Rdata",sep=""))

## visualization
# Choose a set of soft-thresholding powers
powers  =  c(1:30)
#  Call  the  network  topology  analysis  function
sft=pickSoftThreshold(counts.t,powerVector=powers,RsquaredCut=0.9,verbose=5) #power must be between 1 and 30.

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

# create TOM
TOM =TOMsimilarityFromExpr(counts.t,power=sft$powerEstimate) # power=14 shows R^2=0.9
colnames(TOM)=genes
rownames(TOM)=genes
head(TOM)

# extract top hub genes
index.sub=is.element(genes, top.hub)
subTOM=TOM[index.sub,index.sub]

# only strong interaction is shown
h.subTOM = (subTOM>0.1)*subTOM # only > 0.1 TOM will be shown in network

subnet=graph.adjacency(h.subTOM,mode="undirected",weighted=TRUE,diag=FALSE)
summary(subnet)

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
pdf("top.hub.100.pdf", useDingbats=FALSE) 
plot(subnet, layout=layout.graphopt, vertex.size=v.size, vertex.frame.color=NA,vertex.label=v.label, vertex.label.cex=0.05,edge.color="gray57", edge.width=E(subnet)$weight*0.1)
dev.off()
