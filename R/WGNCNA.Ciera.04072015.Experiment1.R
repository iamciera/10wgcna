#WGNCNA.Ciera.04072015.tutorialGuided.R

##Experiment 1

In order to visualize what is going on in the network and understand if co-expression is occuring in other datasets, I will attempt to build adjaceny network and color by SOM cluster.  


```{r}
## library
library(WGCNA)
library (igraph)
library(ggplot2)
library(reshape)

##Setting up the dataset
counts <- read.delim("../data/GSE45774_rpkm_all.txt", header = TRUE)
colnames(counts)[1] <- "Gene_ID"
head(counts)
dim(counts)


#Read in lists of genes from SOM
#SOM_analysis9.5.csv is the large WT only analysis.

SOM <- read.csv("../../08SOM/lcmSOM/data/SOM_analysis9.5forNetwork.csv", header = TRUE) 

#Subset only the columns that specify clusters and genes
SOMsub <- SOM[,c(2,21,22)] 

dim(SOM) # there are only 4,618 genes because I only used the most differentially expressed.

#Subsetted for the interesting genes

head(SOMsub)

#only upregulated in MBR
cluster35 <- subset(SOMsub, som.unit.classif == "35") 

#only upregulated in rachis
cluster17 <- subset(SOMsub, som.unit.classif == "17")

#combine
subclusters <- rbind(cluster35, cluster17)

#isolate only gene names
SOMclusters <- as.data.frame(subclusters[,1])
colnames(SOMclusters)[1] <- "Gene_ID"
dim(SOMclusters) # these are all the genes, there are 144. 

#Now merge with Dan's table to get only genes I am interested in.
#First rename 1st column in counts to gene for merging

merged <- merge(SOMclusters, counts, by = "Gene_ID")

dim(merged) # lost 11 genes that must not have been present in Dan Koenig's data

counts <- merged
counts[is.na(counts)] <- 0
str(counts)

#isolate 
genes <- counts[,1] # for later

# transpose all, but the first column (name)
counts.t <- as.data.frame(t(counts[,-1]))
colnames(counts.t) <- genes #bring the names back
head(counts.t)
dim(counts.t)

datExpr <- counts.t

```

##Making the Adjaceny Matrix

###Picking Soft Threshold

The key here is to pick the "lowest power for which the scale-free topology fit
index curve flattens out upon reaching a high value."Note that at power=5, the curve has an elbow or kink, i.e. for this power the scale free topology fit does not improve after increasing the power. From the the image below, I chose 5.

```{r}
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

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
```

##Building the adjacency matrix
```{r}
softPower = 5;
adjacency = adjacency(datExpr, power = softPower);

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM
head(dissTOM)
# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);

# Set the minimum module size, the larger the min. the larger the clusters.
minModuleSize = 20;

# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")


###########################
#Visulization I - From the package.  Heat Map.
#Based off this tutorial: http://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-05-Visualization.pdf

# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7;

# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;

# Call the plot function
sizeGrWindow(9,9)
TOMplot(plotTOM, geneTree, dynamicColors, main = "Network heatmap plot, all genes")

#Plot specific eigengens
MEs = moduleEigengenes(datExpr, dynamicColors)$eigengenes
# Isolate weight from the clinical traits
weight = as.data.frame(datTraits$weight_g);
names(weight) = "weight"




#Annotation Files
annotation1<- read.delim("../../06diffGeneExp/analysis/data/ITAG2.3_all_Arabidopsis_ITAG_annotations.tsv", header=FALSE)  #Changed to the SGN human readable annotation
colnames(annotation1) <- c("ITAG", "SGN_annotation")
annotation2<- read.delim("../../06diffGeneExp/analysis/data/ITAG2.3_all_Arabidopsis_annotated.tsv")
annotation <- merge(annotation1,annotation2, by = "ITAG")
sub.annotation <- annotation[,c(1,4)] #Choose what you want to name them by

# #Merge with genes then call which column of genes you want.
table.genes <- as.data.frame(genes)
colnames(table.genes) <- "ITAG"
 
head(table.genes)

annotation.merge <- merge(table.genes,sub.annotation, by = "ITAG", all.x = TRUE)

#If NA, replace with ITAG column
annotation.merge$gene.name <- ifelse(is.na(annotation.merge$symbol), 
                                     annotation.merge$ITAG, 
                                     annotation.merge$symbol)

genes <-  annotation.merge$gene.name

#Make Ranking

row.names() <- genes
average <- rowMeans(result)
sd <- apply(result,1,function(d)sd(d))
result.n <- cbind(result,average,sd)
result.n <- as.data.frame(result.n)
result.g <- subset(result.n[,11:12])
qplot(average, sd, data=result.g)  #what should I be seeing here
colnames(result.g) <- c("ave.rank","sd.rank")

result.o <- result.g[order(result.g$ave.rank),]
top.hub <- rownames(result.o[1:400,]) # top  hub genes

