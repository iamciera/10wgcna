#Analysis Key

##WGCNA.bootstrap.Ciera.11.18.14.R
This is setting up the map.  

##WGCNA.bootstrap.Ciera.11.21.14.R


##To Do

1. Try without the all the species data
2. Try with different clusters 
3. Need to define cluster number with colors
4. Try with larger Tomato dataset


##Analysis

`WGNCNA.Ciera.04072015.Experiment1`: Experment one.  In order to visualize what is going on in the network and understand if co-expression is occuring in other datasets, I will attempt to build adjaceny network and color by SOM cluster.  

`WGNCNA.Ciera.04062015.tutorialGuided.R`: I decided to start from the begining, since it is too much trouble to troubleshoot Yasu's script.  This my first attempt guided by tutorial.

`WGCNA.bootstrap.Ciera.03.27.15.c`: **Finally got it working again.**

Data: Yasu's
Clusters: basic SOM, cluster 35.

4.06.15 update: So I came back after a week to figure out what the problem is.  So basically I kept getting this "character" error and I traced it to the transformation `t()` step. I did some work and added another step after the merge to get the `t()` to work again.  Basically, I was forcing the dataframe to transform all wierd, so I need to do the transformation without the gene names and add them later.  I still don't quite understand why it worked before though. 

`WGCNA.bootstrap.Ciera.03.27.15.b`

This is one more attempt to get a Network with some usable information. This time I am going to use Yasus data.  Want to make network with *JUST cluster 35*.

`WGCNA.bootstrap.Ciera.03.27.15`

This is one more attempt to get a Network with some usable information. This time I am going to 
use Dan Koenig's data. Now I am verifying what I did in WGNCNA.bootstrap.Ciera03.26.15.R and also trying to understand why many of the genes are not used in the analysis.

*This will take too Long, maybe run at night.*

`WGCNA.bootstrap.Ciera.03.26.15` #This is one more attempt to get a Network with some usable information. This time I am going to use Dan Koenig's data. This is works to produce a network with three clusters, although it needs to be re-done.  The merging conflict is still present, but works to make a graph.  If you fix the merging conflict you get over 4,000 genes, I have not tried to run this, as it will take a long time!

`WGCNA.bootstrap.Ciera.11.21.14` - Not useful

`WGCNA.bootstrap.Ciera.11.18.14.R`

Approach: The approach is to take the interesting clusters from the SOM and use Yasu's data to inform which genes are 1. major hubs for leaf primordia regulation and 2. Gene clusters of genes that could be working together. 

`WGCNA.bootstrap.Ciera.11.21.14.tutorialGuided`: An attempt to follow the guide, but didn't go so well.


Steps:
1. Input Yasu's data `pnas.1402835111.sd01.csv`
2. Input my data  `SOM_analysis9.5.csv` 
3. Unput my data `superSOM_analysis8.csv`
4. Merge with Yasu's data
5. Perform Network analysis


##Results

The main problem is that after merging I go from 758 to only 152. It seems that those are the only genes that are in common between the genes I find interesting and Yasu's which is just insane. 

I need to look into why.  Maybe use larger data set 




##References

[The Generalized Topological Overlap Matrix For
Detecting Modules in Gene Networks
](http://labs.genetics.ucla.edu/horvath/GTOM/old/GTOM_tech_report.pdf) 