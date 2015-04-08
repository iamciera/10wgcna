#WGNCNA.Ciera.04062015.tutorialGuided.R
#Based on this tutorial: 
#http://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-02-networkConstr-man.pdf

## library
library(WGCNA)
options(stringsAsFactors  =  FALSE)
#enableWGCNAThreads()
ALLOW_WGCNAT_THREADS = 4
library (igraph)
library(ggplot2)
library(reshape)
library(flashClust)

#library(rgl)
#library(tcltk2)

############Setting up the dataset
counts <- read.delim("../data/GSE45774_rpkm_all.txt", header = TRUE)
colnames(counts)[1] <- "Gene_ID"
head(counts)
dim(counts)


#Read in lists of genes from both SOM and superSOM
#SOM_analysis9.5.csv is the large WT only analysis.

SOM <- read.csv("../../08SOM/lcmSOM/data/SOM_analysis9.5forNetwork.csv", header = TRUE) 

#Subset only the columns that specify clusters and genes
SOMsub <- SOM[,c(2,21,22)] 



#Subsetted for the interesting genes

head(SOMsub)

cluster35 <- subset(SOMsub, som.unit.classif == "35")

cluster20 <- subset(SOMsub, som.unit.classif == "20")

cluster5 <- subset(SOMsub, som.unit.classif == "5")

subclusters <- rbind(cluster35, cluster20, cluster5)


#isolate only gene names
SOMclusters <- as.data.frame(subclusters[,1])
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

datExpr <- counts.t

###################################################
#Making the Adjaceny Matrix

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

softPower = 6;
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

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
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
# Add the weight to existing module eigengenes



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

