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
colnames(SOM)
head(SOM)
#Subset only the columns that specify clusters and genes
SOMsub <- SOM[,c(2,21,22)] 
names(SOMsub)


#Subsetted for the interesting genes

head(SOMsub)

cluster35 <- subset(SOMsub, som.unit.classif == "35")
dim(cluster35)

cluster20 <- subset(SOMsub, som.unit.classif == "20")
dim(cluster20)

cluster5 <- subset(SOMsub, som.unit.classif == "5")

subclusters <- rbind(cluster17, cluster20, cluster5)


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
dim(counts.t)


head(counts.t)

#####################

#Guided by the WGCNA tutorial

head(counts.t[,1:4])

gsg = goodSamplesGenes(counts.t, verbose = 3)
gsg$allOK


sampleTree = flashClust(dist(counts.t), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

# Plot a line to show the cut
abline(h = 1000, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 1000, minSize = 10)
table(clust)

# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

###############################
# One-step network construction and module detection

net = blockwiseModules(datExpr, power = 6,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "cieraTOM",
                       verbose = 3)

