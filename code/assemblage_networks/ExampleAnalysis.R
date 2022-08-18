# This example is adapted from the Supplement of Oldham, Horvath, Geschwind
# (2006). The first setwd command should be modified to suit your desired
# working directory.

# Set the working directory of the R session by using the following command:
# setwd("C:/Documents and Settings/plangfelder/My
# Documents/Work/TreeCut/PublicFunction/1.02") Note that we use / instead of \
# in the path. Read in the R libraries:
library(MASS) # standard, no need to install
library(class) # standard, no need to install
library(cluster)
library(sma) # install it for the function plot.mat
library(impute) # install it for imputing missing value
library(Hmisc) # probably you won?t need this
library(splines) # probably, you won?t need this

#Memory
# Check the maximum memory that can be allocated:
# memory.size(TRUE)/1024
# Increase the available memory:
# memory.limit(size=4000)
# Read in the custom network functions:

source("NetworkFunctions-ExampleAnalysis.R");

# Custom libraries:

library(dynamicTreeCut)
library(moduleColor)

## Read in array data:

file = bzfile(description = "Dataset 1 (network construction).csv.bz2")
dat1 = read.csv(file = file, header = T)
attach(dat1)

dim(dat1)
# This data frame contains the gene expression data.
# By our convention, columns are genes and rows are samples.
datExpr = data.frame(t(dat1[Brain_variant_H > 0, 2:39]))
dim(datExpr)
dimnames(datExpr)[[1]]
indexHuman = c(19:36)
indexChimp = c(1:18)

printFlush("Calculating connectivities...")


powerHuman = 9
# Here we use the same power as was used for the human network to facilitate the
# comparison of the two networks. Again, by playing around with different
# choices, one can easily verify that our findings are highly robust with
# respect to the choice of the power.
powerChimp = powerHuman

# Calculation of the network (adjacency matrix) by raising the absolute value of
# the correlation matrix to a power (soft-thresholding with the power adjacency
# function). Human network:
AdjMatHuman = abs(cor(datExpr[indexHuman, ] , use = "p")) ^ powerHuman
diag(AdjMatHuman) = 0
# Chimp network:
AdjMatChimp = abs(cor(datExpr[indexChimp,] , use = "p")) ^ powerChimp
diag(AdjMatChimp) = 0

## Calculation of the whole network connectivity k:
ConnectivityHuman <- apply(AdjMatHuman, 1, sum)
ConnectivityChimp <- apply(AdjMatChimp, 1, sum)
ConnectivityHuman = ConnectivityHuman / max(ConnectivityHuman)
ConnectivityChimp = ConnectivityChimp / max(ConnectivityChimp)

# As a pre-processing step towards module construction, we restrict the network
# to genes with reasonably high connectivity. This does not lead to a big loss
# of information since module genes tend to have high connectivity (7). Toward
# this end, consider the median connectivity in human and chimp:
median(ConnectivityHuman)

# 0.09753368 This motivates us to restrict the analysis to genes with k > 0.1 in
# either human or # chimp:
minconnections = .1
rest1 = ConnectivityChimp > minconnections |
  ConnectivityHuman > minconnections
table(rest1)
AdjMatChimprest = AdjMatChimp[rest1, rest1]
AdjMatHumanrest = AdjMatHuman[rest1, rest1]
rm(AdjMatChimp)
rm(AdjMatHuman)
collect_garbage()

#Module Construction 

# The topological overlap of two nodes reflects their similarity in terms of the
# commonality of the nodes they connect to, see (6, 10).

# Creating distance matrices based upon the topological overlap of 2241 genes
# for humans and chimpanzees: 

# distTOMChimp <- TOMdist1(AdjMatChimprest)

printFlush("Calculating TOM...")
distTOMHuman <- TOMdist1(AdjMatHumanrest)
collect_garbage()

# To group genes with coherent expression profiles into modules, we use average
# linkage hierarchical clustering, which uses the topological overlap measure as
# dissimilarity.

# Performing average linkage hierarchical clustering using these distance
# matrices:
printFlush("Clustering TOM...")
hierTOMHuman <- hclust(as.dist(distTOMHuman), method = "average")
# hierTOMChimp <- hclust(as.dist(distTOMChimp),method="average")

colorh1 = as.character(modulecolor2(hierTOMHuman, h1 = .95, minsize1 = 30))
table(colorh1)

DynamicColor1 = labels2colors(
  cutreeDynamic(
    hierTOMHuman,
    cutHeight = 0.97,
    minClusterSize = 30,
    method = "tree",
    deepSplit = TRUE
  )
)

DynamicColor2 = labels2colors(
  cutreeDynamic(
    hierTOMHuman,
    cutHeight = 0.97,
    minClusterSize = 30,
    method = "tree",
    deepSplit = FALSE
  )
)

ClustColor1 = labels2colors(
  cutreeDynamic(
    dendro = hierTOMHuman,
    minClusterSize = 30,
    cutHeight = 0.97,
    method = "hybrid",
    maxCoreScatter = 0.75,
    minGap = 0.25,
    pamStage = TRUE,
    distM = distTOMHuman,
    useMedoids = FALSE,
    maxDistToLabel = 0.90,
    respectSmallClusters = TRUE
  )
)

Clusters1 = cutreeHybrid(
  dendro = hierTOMHuman,
  minClusterSize = 30,
  cutHeight = 0.97,
  maxCoreScatter = 0.75,
  minGap = 0.250,
  pamStage = TRUE,
  distM = distTOMHuman,
  useMedoids = FALSE,
  maxDistToLabel = 0.90,
  respectSmallClusters = TRUE
)

CoreColor1 = labels2colors(Clusters1$cores);

ClustColor2 = labels2colors(
  cutreeDynamic(
    dendro = hierTOMHuman,
    minClusterSize = 30,
    cutHeight = 0.97,
    method = "hybrid",
    maxCoreScatter = 0.95,
    minGap = 0.050,
    pamStage = TRUE,
    distM = distTOMHuman,
    useMedoids = FALSE,
    maxDistToLabel = 0.90,
    respectSmallClusters = TRUE
  )
)

Clusters2 = cutreeHybrid(
  dendro = hierTOMHuman,
  minClusterSize = 30,
  cutHeight = 0.97,
  maxCoreScatter = 0.95,
  minGap = 0.050,
  pamStage = TRUE,
  distM = distTOMHuman,
  useMedoids = FALSE,
  maxDistToLabel = 0.90,
  respectSmallClusters = TRUE
)

CoreColor2 = labels2colors(Clusters2$cores);

AutoColor = NULL;

for (deepSplit in 0:3)
  AutoColor = cbind(AutoColor, labels2colors(
    cutreeDynamic(
      dendro = hierTOMHuman,
      minClusterSize = 30 - 3 * deepSplit,
      cutHeight = 0.97,
      method = "hybrid",
      deepSplit = deepSplit,
      distM = distTOMHuman
    )
  ))

AutoLabels = paste("Hybrid 'auto': dS =", c(0:3))

par(mfrow = c(2, 1))
par(cex = 1.4)
par(mar = c(0, 8.5, 2, 0))

plot(
  hierTOMHuman,
  labels = F,
  main = "Hierarchical dendrogram and module colors",
  sub = "",
  xlab = ""
)

par(mar = c(1, 8.5, 0, 0))

plotHclustColors(
  hierTOMHuman,
  cbind(colorh1, DynamicColor1, DynamicColor2, #ClustColor1, CoreColor1,
        #ClustColor2, CoreColor2,
        AutoColor),
  c("Static", "Dynamic Tree (dS)", "Dynamic Tree (No dS)", #"Dynamic Hybrid 1",
    #"Cores in DHyb 1", "Dynamic Hybrid 2", "Cores in DHyb 2",
    AutoLabels),
  main = ""
)

