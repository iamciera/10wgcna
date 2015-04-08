## WGCNA with bootstrap
## 2014 1 21
install.packages("WGCNA")
source("http://bioconductor.org/biocLite.R")
biocLite("impute")
biocLite("igraph")
## library
library(WGCNA)
options(stringsAsFactors  =  FALSE)
enableWGCNAThreads()
library (igraph)
library(ggplot2)
library(reshape)

## Read in data: dodder cluster 5 with indicidual biological replicates
## because WGCNA with soft threshold need > 12 samples

# counts <- read.csv("counts.csv",row.names=1)
# dim(counts) #1498   40
# head(counts)

  ann=read.csv("SL_Annotation_20120920.csv",header=T)
  counts=read.csv("IL_normalized_estimated_counts.csv",na.strings="N/A")
  eqtl=read.csv("Clustering network Cluster 24.csv",header=T)
  

##-----##
eqtl.subset=eqtl
	eqtl.subset=eqtl[which(eqtl$colors==4),]

		
    eqtlcounts=counts[which(counts$X %in% eqtl.subset$X.1),]
    rownames(eqtlcounts)=eqtlcounts$X
    eqtlcounts=eqtlcounts[,2:312]
    genes=rownames(eqtlcounts)
    countsadj=matrix(nrow=nrow(eqtlcounts),ncol=ncol(eqtlcounts))

    for(i in 1:nrow(eqtlcounts)){
      for(j in 1:ncol(eqtlcounts)){
        countsadj[i,j]=as.numeric(eqtlcounts[i,j])
      }
      print(i)
    }

# for(i in 1:nrow(countsadj)){
#   for(j in 1:ncol(countsadj)){
# countsadj[i,j]=countsadj[i,j]/countsadj[i,75]
# }
# }

rownames(countsadj)=rownames(eqtlcounts); colnames(countsadj)=colnames(eqtlcounts)
counts.t=t(countsadj)

counts.t=log2(counts.t)

for(i in 1:nrow(counts.t)){
  for(j in 1:ncol(counts.t)){
    if(counts.t[i,j]=="-Inf"){counts.t[i,j]=0}
  }
}

#counts.t=data.Normalization (counts.t,type="n1",normalization="row")
# genes=rownames(counts)
# counts.t <- t(counts)
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
    sft.b = pickSoftThreshold(Y.b,  powerVector=powers, RsquaredCut=0.8, verbose  =  5)
    sft.power = sft.b$powerEstimate
    print(sft.power)
  }
  
  print(sft.power)
  ##TOM
  TOM.b = TOMsimilarityFromExpr(Y.b,power=sft.power) #omega TOM-based connectivity
  hub.b = rowSums(TOM.b)
#   adj.b = adjacency(Y.b,power=sft.b$powerEstimate)
#   hub.b = rowSums(adj.b) #k connectivity
  
  result[,i]<-rank(-hub.b)
}

row.names(result) <- genes
average <- rowMeans(result)
sd <- apply(result,1,function(d)sd(d))
result.n <- cbind(result,average,sd)
result.n <- as.data.frame(result.n)
result.g <- subset(result.n[,c("average","sd")])
qplot(average, sd, data=result.g) 
colnames(result.g) <- c("ave.rank","sd.rank")

# gene annotation
#annotation <- read.csv("Dodder_cdhitest95_TAIR_Blast2GO_annotation.csv", header=TRUE)
#result.a<-merge(result.g,annotation, by.x="row.names", by.y="Sequence_name", all.x=T,sort=F)
#write.csv(result.a,paste("boot",B,".WGCNA.hub.csv",sep=""))

result.o <- result.g[order(result.g$ave.rank),]
top.hub <- rownames(result.o[1:800,]) # top 100 hub genes
write.csv(result.o,file=paste("GOSEQ list.csv"))

rownames(result.o)=substring(rownames(result.o),1,14)
ann.reduced=ann[which(ann$ITAG %in% rownames(result.o)),c("ITAG","Description")]
colnames(ann.reduced)=c("itag","Description")
result.o=merge(result.o,ann.reduced,by.x="row.names",by.y="itag",sort=F)
write.csv(result.o,file=paste("hub list.csv"))

# # save
# save.image(file=paste("boot",B,".WGCNA.Rdata",sep=""))

## visualization
# Choose a set of soft-thresholding powers
powers  =  c(1:30)
#  Call  the  network  topology  analysis  function
sft=pickSoftThreshold(counts.t,powerVector=powers,RsquaredCut=0.75,verbose=5) #power must be between 1 and 30.
print(sft$powerEstimate)

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

between <- fastgreedy.community(subnet,modularity=T,membership=T)

##-----##
community=as.data.frame(matrix(nrow=0,ncol=3))
community=as.matrix(membership(between)); colnames(community)=c("communities")
write.csv(community,file=paste("Community Membership.csv"))
#between.a<-merge(between,annotation, by.x="row.names", by.y="Sequence_name", all.x=T,sort=F)
#write.csv(between.a,"betweenness.csv")
# head(between[order(-between)])
##-----##

names=as.data.frame(rownames(result.o))
names$gene=result.o[,1]; colnames(names)=c("number","gene")
top.hub=as.data.frame(top.hub)

png(file=paste("subcom top.png"),height=1080,width=1920)
plot(subnet,layout=layout.auto,vertex.frame.color="black",vertex.label.cex=0.75,vertex.label.color="dark green",vertex.size=2,vertex.color=membership(between))
dev.off()

png(file=paste("subcom top no label.png"),height=1080,width=1920)
plot(subnet,layout=layout.auto,vertex.frame.color="black",vertex.label=NA,vertex.size=2,vertex.color=membership(between))
dev.off()

##------##

#Junk
##-----##
counts.t.module=community

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
abline(h=0.70,col="red")
#  Mean  connectivity  as  a  function  of  the  soft-thresholding  power
plot(sft$fitIndices[,1],  sft$fitIndices[,5],
     xlab="Soft  Threshold  (power)",ylab="Mean  Connectivity",  type="n",
     main  =  paste("Mean  connectivity"))
text(sft$fitIndices[,1],  sft$fitIndices[,5],  labels=powers,  cex=cex1,col="red")
dev.off()