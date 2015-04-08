#Using WGNCA to build co-expression network

The purpose of performing network analysis is to test if the biologically relevant clusters found in SOM analysis are 1. confirmed with external data 2. use this external data to strengthen and deepen our understanding of the co-expression of these genes. 

The main way I am approaching these aims is to extract the relevant genes from SOM analysis and use the count data from either [Yasu's paper](http://www.pnas.org/content/111/25/E2616.short) and [Dan Koenig's paper](http://www.pnas.org/content/110/28/E2655.short).

There are many ways to approach these tasks, but the main variable to test are:

1. Different Data sets (Dan's and Yasu's)
2. Limit clusters that are relevant (starting gene number)
3. Limit samples within the data set (isolate samples: species, dev. stage, ect.)
4. WGCNA conditions

##Analysis Key

###Analysis

`WGNCNA.Ciera.04072015.Experiment1`: Experiment one.  In order to visualize what is going on in the network and understand if co-expression is occurring in other datasets, I will attempt to build adjacency network and color by SOM cluster.  

`WGNCNA.Ciera.04062015.tutorialGuided.R`: I decided to start from the beginning, since it is too much trouble to troubleshoot Yasu's script.  This my first attempt guided by tutorial.

`WGCNA.bootstrap.Ciera.03.27.15.c`: **Finally got it working again.**

Data: Yasu's
Clusters: basic SOM, cluster 35.

4.06.15 update: So I came back after a week to figure out what the problem is.  So basically I kept getting this "character" error and I traced it to the transformation `t()` step. I did some work and added another step after the merge to get the `t()` to work again.  Basically, I was forcing the data frame to transform all weird, so I need to do the transformation without the gene names and add them later.  I still don't quite understand why it worked before though. 

`WGCNA.bootstrap.Ciera.03.27.15.b`

This is one more attempt to get a Network with some usable information. This time I am going to use Yasu's data.  Want to make network with *JUST cluster 35*.

`WGCNA.bootstrap.Ciera.03.27.15`

This is one more attempt to get a Network with some usable information. This time I am going to 
use Dan Koenig's data. Now I am verifying what I did in WGNCNA.bootstrap.Ciera03.26.15.R and also trying to understand why many of the genes are not used in the analysis.

*This will take too Long, maybe run at night.*

`WGCNA.bootstrap.Ciera.03.26.15` #This is one more attempt to get a Network with some usable information. This time I am going to use Dan Koenig's data. This is works to produce a network with three clusters, although it needs to be re-done.  The merging conflict is still present, but works to make a graph.  If you fix the merging conflict you get over 4,000 genes, I have not tried to run this, as it will take a long time!

`WGCNA.bootstrap.Ciera.11.21.14` - Not useful

`WGCNA.bootstrap.Ciera.11.18.14.R`

Approach: The approach is to take the interesting clusters from the SOM and use Yasu's data to inform which genes are 1. major hubs for leaf primordia regulation and 2. Gene clusters of genes that could be working together. 

`WGCNA.bootstrap.Ciera.11.21.14.tutorialGuided`: An attempt to follow the guide, but didn't go so well.

##References

[The Generalized Topological Overlap Matrix For
Detecting Modules in Gene Networks
](http://labs.genetics.ucla.edu/horvath/GTOM/old/GTOM_tech_report.pdf) 