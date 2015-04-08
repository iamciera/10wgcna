## WGCNA with bootstrap
## 3.27.16.d
## Ciera Martinez

#This is one more attempt to get a Network with some usable information. This time I am going to 
#use Dan Koenig's data.  want to make network with JUST cluster 35

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


## Yasu 
counts <- read.csv("../data/pnas.1402835111.sd01.csv")
colnames(counts)[1] <- "Gene_ID"
head(counts)
dim(counts)

#at this point do I need to get rid of all the samples I don't need
#First I need to subset based on the genes I am interested in form my analysis. 

#Read in lists of genes from both SOM and superSOM
#SOM_analysis9.5.csv is the large WT only analysis.
#The problem is here, I was unsing the analysis9.5, which is truncated to fit with Yasu's curated data. 
SOM <- read.csv("../../08SOM/lcmSOM/data/SOM_analysis9.5forNetwork.csv", header = TRUE) 
colnames(SOM)
head(SOM)
#Subset only the columns that specify clusters and genes
SOMsub <- SOM[,c(2,21,22)] 
names(SOMsub)


#OMG I never subseted for the interesting genes?
# Maybe I icluded everything to make sure I had enough awesome genes?  I want to run with just 
#cluster 35.  Whoop Whoop.
head(SOMsub)

cluster35 <- subset(SOMsub, som.unit.classif == "35")
dim(cluster35)

cluster17 <- subset(SOMsub, som.unit.classif == "17")
dim(cluster17)


#isolate only gene names
SOMclusters <- as.data.frame(cluster35[,1])
colnames(SOMclusters)[1] <- "Gene_ID"
head(SOMclusters)

#remove duplicates
dim(SOMclusters)
SOMclusters <- unique(SOMclusters)
dim(SOMclusters)

#Now merge with Dan's table to get only genes I am interested in.
#First rename 1st column in counts to gene for merging

#Merge
dim(counts)
head(counts)
dim(SOMclusters)

merged <- merge(SOMclusters, counts, by = "Gene_ID")


#set row and column names R
#remove rows that the have duplicate gene names.
dim(merged)
countsUniq <- unique(merged)
dim(countsUniq)
#no duplicates


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

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
sft = pickSoftThreshold(counts.t, powerVector = powers, verbose = 5)

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")



#in Yasu's script, why is this used?
#counts[, c(2:64)] <- sapply(counts[, c(2:64)], as.numeric)
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
top.hub <- rownames(result.o[1:200,]) # top  hub genes

# save
save.image(file=paste("boot",B,".WGCNA.Rdata",sep=""))

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
pdf("top.hub.test.032715.a.pdf", useDingbats=FALSE) 
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

