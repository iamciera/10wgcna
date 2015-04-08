# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored.
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = ".";
setwd(workingDir);
# Load WGCNA package
library(WGCNA)
library(cluster)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);


# Here are input parameters of the simulation model
# number of samples or microarrays in the training data
no.obs=50
# now we specify the true measures of eigengene significance
# recall that ESturquoise=cor(y,MEturquoise)
ESturquoise=0; ESbrown= -.6;
ESgreen=.6;ESyellow=0
# Note that we dont specify the eigengene significance of the blue module
# since it is highly correlated with the turquoise module.
ESvector=c(ESturquoise,ESbrown,ESgreen,ESyellow)
# number of genes
nGenes1=3000
# proportion of genes in the turquoise, blue, brown, green, and yellow module #respectively.
simulateProportions1=c(0.2,0.15, 0.08, 0.06, 0.04)
# Note that the proportions dont add up to 1. The remaining genes will be colored grey,
# ie the grey genes are non-module genes.
# set the seed of the random number generator. As a homework exercise change this seed.
set.seed(1)
#Step 1: simulate a module eigengene network.
# Training Data Set I
MEgreen=rnorm(no.obs)
scaledy=MEgreen*ESgreen+sqrt(1-ESgreen^2)*rnorm(no.obs)
y=ifelse( scaledy>median(scaledy),2,1)
MEturquoise= ESturquoise*scaledy+sqrt(1-ESturquoise^2)*rnorm(no.obs)
# we simulate a strong dependence between MEblue and MEturquoise
MEblue= .6*MEturquoise+ sqrt(1-.6^2) *rnorm(no.obs)
MEbrown= ESbrown*scaledy+sqrt(1-ESbrown^2)*rnorm(no.obs)
MEyellow= ESyellow*scaledy+sqrt(1-ESyellow^2)*rnorm(no.obs)
ModuleEigengeneNetwork1=data.frame(y,MEturquoise,MEblue,MEbrown,MEgreen, MEyellow)

dat1=simulateDatExpr5Modules(MEturquoise=ModuleEigengeneNetwork1$MEturquoise,
                             MEblue=ModuleEigengeneNetwork1$MEblue,
                             MEbrown=ModuleEigengeneNetwork1$MEbrown,
                             MEyellow=ModuleEigengeneNetwork1$MEyellow,
                             MEgreen=ModuleEigengeneNetwork1$MEgreen,
                             nGenes=nGenes1,
                             simulateProportions=simulateProportions1)

datExpr = dat1$datExpr;
truemodule = dat1$truemodule;
datME = dat1$datME;
attach(ModuleEigengeneNetwork1)

table(truemodule)
dim(datExpr)

datExpr=data.frame(datExpr)
ArrayName=paste("Sample",1:dim(datExpr)[[1]], sep="" )
# The following code is useful for outputting the simulated data
GeneName=paste("Gene",1:dim(datExpr)[[2]], sep="" )
dimnames(datExpr)[[1]]=ArrayName
dimnames(datExpr)[[2]]=GeneName

rm(dat1); collectGarbage();
# The following command will save all variables defined in the current session.
save.image("Simulated-dataSimulation.RData");

rm(dat1); collectGarbage();
# The following command will save all variables defined in the current session.
save.image("Simulated-dataSimulation.RData");


# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored.
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = ".";
setwd(workingDir);
# Load WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the previously saved data
load("./")
load("Simulated-dataSimulation.RData");
attach(ModuleEigengeneNetwork1)

GS1= as.numeric(cor(y, datExpr, use="p"))
head(datExpr)
# Network terminology: GS1 will be referred to as signed gene significance measure
p.Standard=corPvalueFisher(GS1, nSamples =length(y) )
# since the q-value function has problems with missing data, we use the following trick
p.Standard2=p.Standard
p.Standard2[is.na(p.Standard)]=1
q.Standard=qvalue(p.Standard2)$qvalues
# Form a data frame to hold the results
StandardGeneScreeningResults=data.frame(GeneName,PearsonCorrelation=GS1, p.Standard, q.Standard)

head(StandardGeneScreeningResults)

NoiseGeneIndicator=is.element( truemodule, c("turquoise", "blue", "yellow", "grey"))+.0
SignalGeneIndicator=1-NoiseGeneIndicator

table(q.Standard<.20)
mean(NoiseGeneIndicator[q.Standard<=0.20])

save.image(file = "Simulated-StandardScreening.RData")

######################################
#5. Construction of a weighted gene co-expression network and network modules

# here we define the adjacency matrix using soft thresholding with beta=6
ADJ1=abs(cor(datExpr,use="p"))^6
dim(datExpr)

# When you have relatively few genes (<5000) use the following code
k=as.vector(apply(ADJ1,2,sum, na.rm=T))

# Plot a histogram of k and a scale free topology plot
sizeGrWindow(10,5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k, main="Check scale free topology\n")

#5.c Comparing various module detection methods
# Turn adjacency into a measure of dissimilarity
dissADJ=1-ADJ1

dissTOM=TOMdist(ADJ1)
collectGarbage()

#5.c.3 Partitioning Around Medoids based on adjacency

pam4=pam(as.dist(dissADJ), 4)
pam5=pam(as.dist(dissADJ), 5)
pam6=pam(as.dist(dissADJ), 6)
# Cross-tabulte the detected and the true (simulated) module membership:
table(pam4$clustering, truemodule)
table(pam5$clustering, truemodule)
table(pam6$clustering, truemodule)

#5.c.5 Average linkage hierachical clustering with adjacency-based dissimilarity

hierADJ=hclust(as.dist(dissADJ), method="average" )
# Plot the resulting clustering tree together with the true color assignment
sizeGrWindow(10,5);
plotDendroAndColors(hierADJ, colors = data.frame(truemodule), dendroLabels = FALSE, hang = 0.03,
                    main = "Gene hierarchical clustering dendrogram and simulated module colors" )

colorStaticADJ=as.character(cutreeStaticColor(hierADJ, cutHeight=.99, minSize=20))
# Plot the dendrogram with module colors
sizeGrWindow(10,5);
plotDendroAndColors(hierADJ, colors = data.frame(truemodule, colorStaticADJ),
                    dendroLabels = FALSE, abHeight = 0.99,
                    main = "Gene dendrogram and module colors")

#5.c.7 Module definition via dynamic branch cutting methods

branch.number=cutreeDynamic(hierADJ,method="tree")
# This function transforms the branch numbers into colors
colorDynamicADJ=labels2colors(branch.number )

colorDynamicHybridADJ=labels2colors(cutreeDynamic(hierADJ,distM= dissADJ,
                                                  cutHeight = 0.998, deepSplit=2, pamRespectsDendro = FALSE))

# Plot results of all module detection methods together:
sizeGrWindow(10,5)
plotDendroAndColors(dendro = hierADJ,
                    colors=data.frame(truemodule, colorStaticADJ,
                                      colorDynamicADJ, colorDynamicHybridADJ),
                    dendroLabels = FALSE, marAll = c(0.2, 8, 2.7, 0.2),
                    main = "Gene dendrogram and module colors")

##5.c.8 Module definition using the topological overlap based dissimilarity

# Calculate the dendrogram
hierTOM = hclust(as.dist(dissTOM),method="average");
# The reader should vary the height cut-off parameter h1
# (related to the y-axis of dendrogram) in the following
colorStaticTOM = as.character(cutreeStaticColor(hierTOM, cutHeight=.99, minSize=20))
colorDynamicTOM = labels2colors (cutreeDynamic(hierTOM,method="tree"))
colorDynamicHybridTOM = labels2colors(cutreeDynamic(hierTOM, distM= dissTOM , cutHeight = 0.998,
                                                    deepSplit=2, pamRespectsDendro = FALSE))
# Now we plot the results
sizeGrWindow(10,5)
plotDendroAndColors(hierTOM,
                    colors=data.frame(truemodule, colorStaticTOM,
                                      colorDynamicTOM, colorDynamicHybridTOM),
                    dendroLabels = FALSE, marAll = c(1, 8, 3, 1),
                    main = "Gene dendrogram and module colors, TOM dissimilarity")

tabStaticADJ=table(colorStaticADJ,truemodule)
tabStaticTOM=table(colorStaticTOM,truemodule)
tabDynamicADJ=table(colorDynamicADJ, truemodule)
tabDynamicTOM=table(colorDynamicTOM,truemodule)
tabDynamicHybridADJ =table(colorDynamicHybridADJ,truemodule)
tabDynamicHybridTOM =table(colorDynamicHybridTOM,truemodule)
randIndex(tabStaticADJ,adjust=F)
randIndex(tabStaticTOM,adjust=F)
randIndex(tabDynamicADJ,adjust=F)
randIndex(tabDynamicTOM,adjust=F)
randIndex(tabDynamicHybridADJ ,adjust=F)
randIndex(tabDynamicHybridTOM ,adjust=F)

colorh1= colorDynamicHybridTOM
# remove the dissimilarities, adjacency matrices etc to free up space
rm(ADJ1); rm(dissADJ);
collectGarbage()
save.image("Simulated-NetworkConstruction.RData")

# Load the previously saved data
load("Simulated-RelatingToExt.RData");
