# PL: removed cutreeDynamic as this is now implemented in library treecut.

# This document contains functions that we find useful for gene co-expression network analysis
# We recommend to read the tutorial first.
# The functions were written by Steve Horvath, Jun Dong, Peter Langfelder, Andy Yip, Bin Zhang
# To cite this code or the statistical methods please use the following 2 references:

#1) Horvath S, Zhang B, Carlson M, Lu KV, Zhu S, Felciano RM, Laurance MF, Zhao W, Shu, Q, Lee Y, 
# Scheck AC, Liau LM, Wu H, Geschwind DH, Febbo PG, Kornblum HI, Cloughesy TF, Nelson SF, Mischel PS (2006) 
#"Analysis of Oncogenic Signaling Networks in Glioblastoma Identifies ASPM as a Novel Molecular Target", 
# PNAS | November 14, 2006 | vol. 103 | no. 46 | 17402-17407

#2) Zhang B, Horvath S (2005) "A General Framework for Weighted Gene Co-Expression Network Analysis", 
#Statistical Applications in Genetics and Molecular Biology: Vol. 4: No. 1, Article 17 


# Technical Report and software code at: www.genetics.ucla.edu/labs/horvath/CoexpressionNetwork.

# CONTENTS  
# This document contains function for carrying out the following tasks
# A) Assessing scale free topology and choosing the parameters of the adjacency function
#    using the scale free topology criterion (Zhang and Horvath 05)
# B) Computing the topological overlap matrix 
# C) Defining gene modules using clustering procedures
# D) Summing up modules by their first principal component (first eigengene)
# E) Relating a measure of gene significance to the modules 
# F) Carrying out a within module analysis (computing intramodular connectivity) 
#    and relating intramodular connectivity to gene significance.
# G) Miscellaneous other functions, e.g. for computing the cluster coefficient.
# H) Network functions from Bin Zhang and Peter Langfelder for dynamic tree cutting 
     #of a hierarchical clustering tree (dendrogram)
# I) Peter Langfelder's functions for merging modules, e.g. when their module eigengenes have a high correlation
# J) General statistical functions, e.g. for scatterplots

# The functions below are organized according into categories A-J.

#####################################################################################################
################################################################################################################################
# A) Assessing scale free topology and choosing the parameters of the adjacency function
#    using the scale free topology criterion (Zhang and Horvath 05)
#################################################################################################################


# ===================================================
#For hard thresholding, we use the signum (step) function
if(exists("signum1") ) rm(signum1); 
signum1=function(corhelp,tau1)  {
adjmat1= as.matrix(abs(corhelp)>=tau1)
dimnames(adjmat1) <- dimnames(corhelp)
diag(adjmat1) <- 0
adjmat1}

# ===================================================
# For soft thresholding, one can use the sigmoid function 
# But we have focused on the power adjacency function in the tutorial...
if (exists("sigmoid1") ) rm(sigmoid1); sigmoid1=function(ss,mu1=0.8,alpha1=20) {
1/(1+exp(-alpha1*(ss-mu1)))}





#This function is useful for speeding up the connectivity calculation.
#The idea is to partition the adjacency matrix into consecutive baches of a #given size.
#In principle, the larger the batch size the faster is the calculation. But #smaller batchsizes require #less memory...
# Input: gene expression data set where *rows* correspond to microarray samples #and columns correspond to genes. 
# If fewer than MinimumNoSamples contain gene expression information for a given
# gene, then its connectivity is set to missing. 
if(exists("SoftConnectivity")) rm(SoftConnectivity);
SoftConnectivity=function(datE, power=6,batchsize=1500,MinimumNoSamples=10) {
no.genes=dim(datE)[[2]]
no.samples=dim(datE)[[1]]
if (no.genes<no.samples | no.genes<10 | no.samples<5 ) {print("Error: Something seems to be wrong. Make sure that the input data frame has genes as rows and array samples as columns. Alternatively, there could be fewer than 10 genes or fewer than 5 samples. ") } else {
sum1=function(x) sum(x,na.rm=T)
k=rep(NA,no.genes)
no.batches=as.integer(no.genes/ batchsize)
if (no.batches>0) {
for (i in 1:no.batches) {
print(paste("batch number = ", i))
index1=c(1:batchsize)+(i-1)* batchsize
ad1=abs(cor(datE[,index1], datE,use="p"))^power
ad1[is.na(ad1)]=0
k[index1]=apply(ad1,1,sum1)
# If fewer than MinimumNoSamples contain gene expression information for a given
# gene, then we set its connectivity to 0.
NoSamplesAvailable=apply(!is.na(datE[,index1]),2,sum)
k[index1][NoSamplesAvailable< MinimumNoSamples]=NA
} # end of for (i in 1:no.batches
} # end of if (no.batches>0)...
if (no.genes-no.batches*batchsize>0 ) {
restindex=c((no.batches*batchsize+1):no.genes)
ad1=abs(cor(datE[,restindex], datE,use="p"))^power
ad1[is.na(ad1)]=0
k[restindex]=apply(ad1,1,sum1)
NoSamplesAvailable=apply(!is.na(datE[,restindex]),2,sum)
k[restindex][NoSamplesAvailable< MinimumNoSamples]=NA
} # end of if
} # end of else statement
k
} # end of function




# ===================================================
# The function PickHardThreshold can help one to estimate the cut-off value 
# when using the signum (step) function.
# The first column lists the threshold ("cut"),
# the second column lists the corresponding p-value based on the Fisher Transform 
# of the correlation. 
# The third column reports the resulting scale free topology fitting index R^2.
# The fourth column reports the slope of the fitting line, it shoud be negative for 
# biologically meaningul networks.
# The fifth column reports the fitting index for the truncated exponential model. 
# Usually we ignore it.
# The remaining columns list the mean, median and maximum resulting connectivity.
# To pick a hard threshold (cut) with the scale free topology criterion:
# aim for high scale free R^2 (column 3), high connectivity (col 6) and negative slope 
# (around -1, col 4).
# The output is a list with 2 components. The first component lists a sugggested cut-off
# while the second component contains the whole table.
# The removeFirst option removes the first point (k=0, P(k=0)) from the regression fit.
# no.breaks specifies how many intervals used to estimate the frequency p(k) i.e. the no. of points in the 
# scale free topology plot.
if (exists("PickHardThreshold")) rm(PickHardThreshold);
PickHardThreshold=function(datExpr1,RsquaredCut=0.85, cutvector=seq(0.1,0.9,by=0.05) ,removeFirst=FALSE,no.breaks=10) {
no.genes   <- dim(datExpr1)[[2]]
no.genes <- dim(datExpr1)[[2]]
no.samples= dim(datExpr1)[[1]]
colname1=c("Cut","p-value", "scale law R^2", "slope="  ,"truncated R^2","mean(k)","median(k)","max(k)")
datout=data.frame(matrix(666,nrow=length(cutvector),ncol=length(colname1) ))
names(datout)=colname1
datout[,1]=cutvector
for (i in c(1:length(cutvector) ) ){
cut1=cutvector[i]
datout[i,2]=2*(1-pt(sqrt(no.samples-1)*cut1/sqrt(1-cut1^2),no.samples-1))}
if(exists("fun1")) rm(fun1)
fun1=function(x) {
corx=abs(cor(x,datExpr1,use="p"))
out1=rep(NA, length(cutvector) )
for (j in c(1:length(cutvector))) {out1[j]=sum(corx>cutvector[j])}
out1
} # end of fun1
datk=t(apply(datExpr1,2,fun1))
for (i in c(1:length(cutvector) ) ){
nolinkshelp <- datk[,i]-1
cut2=cut(nolinkshelp,no.breaks)
binned.k=tapply(nolinkshelp,cut2,mean)
freq1=as.vector(tapply(nolinkshelp,cut2,length)/length(nolinkshelp))
# The following code corrects for missing values etc
breaks1=seq(from=min(nolinkshelp),to=max(nolinkshelp),length=no.breaks+1)
hist1=hist(nolinkshelp,breaks=breaks1,equidist=F,plot=FALSE,right=TRUE)
binned.k2=hist1$mids
binned.k=ifelse(is.na(binned.k),binned.k2,binned.k)
binned.k=ifelse(binned.k==0,binned.k2,binned.k)
freq1=ifelse(is.na(freq1),0,freq1)
xx= as.vector(log10(binned.k))
if(removeFirst) {freq1=freq1[-1]; xx=xx[-1]}
plot(xx,log10(freq1+.000000001),xlab="log10(k)",ylab="log10(p(k))" )
lm1= lm(as.numeric(log10(freq1+.000000001))~ xx )
lm2=lm(as.numeric(log10(freq1+.000000001))~ xx+I(10^xx) )
datout[i,3]=summary(lm1)$adj.r.squared 
datout[i,4]=summary(lm1)$coefficients[2,1]  
datout[i,5]=summary(lm2)$adj.r.squared
datout[i,6]=mean(nolinkshelp)
datout[i,7]= median(nolinkshelp)
datout[i,8]= max(nolinkshelp) 
} 
datout=signif(datout,3) 
print(data.frame(datout));
# the cut-off is chosen as smallest cut with R^2>RsquaredCut 
ind1=datout[,3]>RsquaredCut
indcut=NA
indcut=ifelse(sum(ind1)>0,min(c(1:length(ind1))[ind1]),indcut)
# this estimates the cut-off value that should be used. 
# Don't trust it. You need to consider slope and mean connectivity as well!
cut.estimate=cutvector[indcut][[1]]
list(cut.estimate, data.frame(datout));
} # end of function











# ===========================================================
# The function PickSoftThreshold allows one to estimate the power parameter when using
# a soft thresholding approach with the use of the power function AF(s)=s^Power
# The function PickSoftThreshold allows one to estimate the power parameter when using
# a soft thresholding approach with the use of the power function AF(s)=s^Power
# The removeFirst option removes the first point (k=1, P(k=1)) from the regression fit.
if (exists("PickSoftThreshold")) rm(PickSoftThreshold);
PickSoftThreshold=function(datExpr1,RsquaredCut=0.85, powervector=c(seq(1,10,by=1),seq(12,20,by=2)),
removeFirst=FALSE,no.breaks=10) {
no.genes   <- dim(datExpr1)[[2]]
no.genes <- dim(datExpr1)[[2]]
no.samples= dim(datExpr1)[[1]]
colname1=c("Power", "scale law R^2" ,"slope", "truncated R^2","mean(k)","median(k)","max(k)")
datout=data.frame(matrix(666,nrow=length(powervector),ncol=length(colname1) ))
names(datout)=colname1
datout[,1]=powervector
if(exists("fun1")) rm(fun1)
fun1=function(x) {
corx=abs(cor(x,datExpr1,use="p"))
out1=rep(NA, length(powervector) )
for (j in c(1:length(powervector))) {out1[j]=sum(corx^powervector[j])}
out1
} # end of fun1
datk=t(apply(datExpr1,2,fun1))
for (i in c(1:length(powervector) ) ){
nolinkshelp <- datk[,i]-1
cut2=cut(nolinkshelp,no.breaks)
binned.k=tapply(nolinkshelp,cut2,mean)
freq1=as.vector(tapply(nolinkshelp,cut2,length)/length(nolinkshelp))
# The following code corrects for missing values etc
breaks1=seq(from=min(nolinkshelp),to=max(nolinkshelp),length=no.breaks+1)
hist1=hist(nolinkshelp,breaks=breaks1,equidist=F,plot=FALSE,right=TRUE)
binned.k2=hist1$mids
binned.k=ifelse(is.na(binned.k),binned.k2,binned.k)
binned.k=ifelse(binned.k==0,binned.k2,binned.k)
freq1=ifelse(is.na(freq1),0,freq1)

xx= as.vector(log10(binned.k))
if(removeFirst) {freq1=freq1[-1]; xx=xx[-1]}
plot(xx,log10(freq1+.000000001),xlab="log10(k)",ylab="log10(p(k))" )
lm1= lm(as.numeric(log10(freq1+.000000001))~ xx )
lm2=lm(as.numeric(log10(freq1+.000000001))~ xx+I(10^xx) )
datout[i,2]=summary(lm1)$adj.r.squared 
datout[i,3]=summary(lm1)$coefficients[2,1]  
datout[i,4]=summary(lm2)$adj.r.squared
datout[i,5]=mean(nolinkshelp)
datout[i,6]= median(nolinkshelp)
datout[i,7]= max(nolinkshelp) 
} 
datout=signif(datout,3) 
print(data.frame(datout));
# the cut-off is chosen as smallest cut with R^2>RsquaredCut 
ind1=datout[,2]>RsquaredCut
indcut=NA
indcut=ifelse(sum(ind1)>0,min(c(1:length(ind1))[ind1]),indcut)
# this estimates the power value that should be used. 
# Don't trust it. You need to consider slope and mean connectivity as well!
power.estimate=powervector[indcut][[1]]
list(power.estimate, data.frame(datout));
}






# ===================================================
# The function ScaleFreePlot1 creates a plot for checking scale free topology
# when truncated1=T is specificed, it provides the R^2 measures for the following
# degree distributions: a) scale free topology, b) log-log R^2 and c) truncated exponential R^2

# The function ScaleFreePlot1 creates a plot for checking scale free topology
if(exists("ScaleFreePlot1")) rm(ScaleFreePlot1) ; 
ScaleFreePlot1=function(kk,no.breaks=10,AF1="" ,truncated1=FALSE, removeFirst=FALSE,cex.lab1=1){
cut1=cut(kk,no.breaks)
binned.k=tapply(kk,cut1,mean)
freq1=tapply(kk,cut1,length)/length(kk)
# The following code corrects for missing values etc
breaks1=seq(from=min(kk),to=max(kk),length=no.breaks+1)
hist1=hist(kk,breaks=breaks1,equidist=F,plot=FALSE,right=TRUE)
binned.k2=hist1$mids
binned.k=ifelse(is.na(binned.k),binned.k2,binned.k)
binned.k=ifelse(binned.k==0,binned.k2,binned.k)
freq1=ifelse(is.na(freq1),0,freq1)
plot(log10(binned.k),log10(freq1+.000000001),xlab="log10(k)",ylab="log10(p(k))",cex.lab=cex.lab1 )
xx= as.vector(log10(binned.k))
if(removeFirst) {freq1=freq1[-1]; xx=xx[-1]}
lm1=lm(as.numeric(log10(freq1+.000000001))~ xx )
lines(xx,predict(lm1),col=1)
OUTPUT=data.frame(ScaleFreeRsquared=round(summary(lm1)$adj.r.squared,2),Slope=round(lm1$coefficients[[2]],2))
if (truncated1==TRUE) { 
lm2=lm(as.numeric(log10(freq1+.000000001))~ xx+I(10^xx) );
OUTPUT=data.frame(ScaleFreeRsquared=round(summary(lm1)$adj.r.squared,2),Slope=round(lm1$coefficients[[2]],2),
TruncatedRsquared=round(summary(lm2)$adj.r.squared,2))
print("the red line corresponds to the truncated exponential fit")
lines(xx,predict(lm2),col=2);
title(paste(AF1, 
", scale free R^2=",as.character(round(summary(lm1)$adj.r.squared,2)),
", slope=", round(lm1$coefficients[[2]],2),
", trunc.R^2=",as.character(round(summary(lm2)$adj.r.squared,2))))} else { 
title(paste(AF1, ", scale R^2=",as.character(round(summary(lm1)$adj.r.squared,2)) , 
", slope=", round(lm1$coefficients[[2]],2)))
}
OUTPUT
} # end of function 









#################################################################################################################
################################################################################################################################
# B) Computing the topological overlap matrix 
#################################################################################################################
#################################################################################################################



# ===================================================
#The function TOMdist1 computes a dissimilarity 
# based on the topological overlap matrix (Ravasz et al)
# Input: an Adjacency matrix with entries in [0,1]
if(exists("TOMdist1")) rm(TOMdist1);
TOMdist1=function(adjmat1, maxADJ=FALSE) {
diag(adjmat1)=0;
adjmat1[is.na(adjmat1)]=0;
maxh1=max(as.dist(adjmat1) ); minh1=min(as.dist(adjmat1) ); 
if (maxh1>1 | minh1 < 0 ) {print(paste("ERROR: the adjacency matrix contains entries that are larger than 1 or smaller than 0!!!, max=",maxh1,", min=",minh1)) } else { 
if (  max(c(as.dist(abs(adjmat1-t(adjmat1)))))>10^(-12)   ) {print("ERROR: non-symmetric adjacency matrix!!!") } else { 
adjmat1= (adjmat1+ t(adjmat1) )/2
kk=apply(adjmat1,2,sum)
maxADJconst=1
if (maxADJ==TRUE) maxADJconst=max(c(as.dist(adjmat1 ))) 
Dhelp1=matrix(kk,ncol=length(kk),nrow=length(kk))
denomTOM= pmin(as.dist(Dhelp1),as.dist(t(Dhelp1)))   +as.dist(maxADJconst-adjmat1); 
gc();gc();
numTOM=as.dist(adjmat1 %*% adjmat1 +adjmat1);
#TOMmatrix=numTOM/denomTOM
# this turns the TOM matrix into a dissimilarity 
out1=1-as.matrix(numTOM/denomTOM) 
diag(out1)=1 
# setting the diagonal to 1 is unconventional (it should be 0)
# but it leads to nicer looking TOM plots... 
out1
}}}





# ===================================================
# This function computes a TOMk dissimilarity
# which generalizes the topological overlap matrix (Ravasz et al)
# Input: an Adjacency matrix with entries in [0,1]
# WARNING:  ONLY FOR UNWEIGHTED NETWORKS, i.e. the adjacency matrix contains binary entries...
# This function is explained in Yip and Horvath (2005)
# http://www.genetics.ucla.edu/labs/horvath/GTOM/
if(exists("TOMkdist1")) rm(TOMkdist1);
TOMkdist1 = function(adjmat1,k=1){
    maxh1=max(as.dist(adjmat1) ); minh1=min(as.dist(adjmat1) );
    if (k!=round(abs(k))) {
        stop("k must be a positive integer!!!", call.=TRUE);}
    if (maxh1>1 | minh1 < 0 ){
        print(paste("ERROR: entries of the adjacency matrix must be between inclusively 0 and 1!!!, max=",maxh1,", min=",minh1))}
    else {
if (  max(c(as.dist(abs(adjmat1-t(adjmat1)))))>0   ) {print("ERROR: non-symmetric adjacency matrix!!!") } else { 

        B <- adjmat1;
        if (k>=2) {
            for (i in 2:k) {
                diag(B) <- diag(B) + 1;
                B = B %*% adjmat1;}}   # this gives the number of paths with length at most k connecting a pair
        B <- (B>0);   # this gives the k-step reachability from a node to another
        diag(B) <- 0;   # exclude each node being its own neighbor
        B <- B %*% B   # this gives the number of common k-step-neighbor that a pair of nodes share

        Nk <- diag(B);
        B <- B +adjmat1;   # numerator
        diag(B) <- 1;
        denomTOM=outer(Nk,Nk,FUN="pmin")+1-adjmat1;
        diag(denomTOM) <- 1;
        1 - B/denomTOM   # this turns the TOM matrix into a dissimilarity
}}
}


# IGNORE THIS function...
# The function TOMdistROW computes the TOM distance of a gene (node)
# with that of all other genes in the network.
# WhichRow is an integer that specifies which row of the adjacency matrix
# corresponds to the gene of interest.
# Output=vector of TOM distances.
if (exists("TOMdistROW") ) rm(TOMdistROW) 
TOMdistROW=function(WhichRow=1, adjmat1, maxADJ=FALSE) {
diag(adjmat1)=0;
maxh1=max(as.dist(adjmat1) ); minh1=min(as.dist(adjmat1) ); 
if (maxh1>1 | minh1 < 0 ) {print(paste("ERROR: the adjacency matrix contains entries that are larger than 1 or smaller than 0!!!, max=",maxh1,", min=",minh1)) } else { 
kk=apply(adjmat1,2,sum)
numTOM=adjmat1[WhichRow,] %*% adjmat1 +adjmat1[WhichRow,]; 
numTOM[WhichRow]=1
maxADJconst=1
if (maxADJ==TRUE) maxADJconst=max(c(as.dist(adjmat1 ))) 
denomTOM=pmin(kk[WhichRow],kk)+maxADJconst-adjmat1[WhichRow,]; denomTOM[WhichRow]=1
#TOMmatrix=numTOM/denomTOM
# this turns the TOM matrix into a dissimilarity 
1-numTOM/denomTOM 
}
}


#####################################################################################################
################################################################################################################################
# C) Defining gene modules using clustering procedures
#####################################################################################################
################################################################################################################################


# ===================================================
#The function modulecolor2 function assigns colors to the observations 
# in the branches of a dendrogram. Now enhanced to handle all possible R colors. 
# Needs GlobalStandardColors.

# we use it to define modules....

if (exists("modulecolor2")) rm(modulecolor2);
modulecolor2=function(hier1, h1=0.9,minsize1=50) {
# here we define modules by using a height cut-off for the branches
labelpred= cutree(hier1,h=h1)
sort1=-sort(-table(labelpred))
modulename= as.numeric(names(sort1))
modulebranch= sort1>minsize1
no.modules=sum(modulebranch)
# now we assume that there are fewer than a certain number of colors
colorcode=GlobalStandardColors
# "grey" means not in any module;
colorhelp=rep("grey",length(labelpred))
if ( no.modules==0 | no.modules >length(colorcode)){ print(paste("The number of modules is problematic. Number of modules = ", as.character(no.modules)))} else { for (i in c(1:no.modules)) {colorhelp=ifelse(labelpred==modulename[i],colorcode[i],colorhelp)};
colorhelp=factor(colorhelp,levels=c(colorcode[1:no.modules],"grey"))
}
factor(colorhelp, levels=unique(colorhelp[hier1$order] ))
}





#------------------------------------------------------------------------------------------------
# Create a barplot of colors given in the couleur parameter ordered according
# to the order of the hierarchical tree hier1. 
# Parameters:
#    hier1:		Hierarchical tree (such as returned by hclust).
#    couleur:		Either a vector of colors for each element in the tree, or a matrix in which every row 
#  			is a set of colors for the corresponding object in the tree. If a matrix, the colors
#			will be plotted starting at the	bottom. 
#    RowLabels:		If given, each color row will be labeled by the corresponding entry in RowLabels; 
#			otherwise the number of the row will be used as a label.
#    cex.RowLabels:	Scale factor for the font size used for labeling rows.
#
# Return value:	None.
#    
# hclustplot1.old implements the case where couleur is a vector, hclustplotn
# the case of a matrix; hclustplot1 is a common wrapper.






#------------------------------------------------------------------------------------------------
#
# hclustplot1, hclustplot1.old, hclustplotn
# The funtion hclusplotn was created by Peter Langfelder.
#------------------------------------------------------------------------------------------------
# Create a barplot of colors given in the couleur parameter ordered according
# to the order of the hierarchical tree hier1. 
# Parameters:
#    hier1:		Hierarchical tree (such as returned by hclust).
#    couleur:		Either a vector of colors for each element in the tree, or a matrix in which every row 
#  			is a set of colors for the corresponding object in the tree. If a matrix, the colors
#			will be plotted starting at the	bottom. 
#    RowLabels:		If given, each color row will be labeled by the corresponding entry in RowLabels; 
#			otherwise the number of the row will be used as a label.
#    cex.RowLabels:	Scale factor for the font size used for labeling rows.
#
# Return value:	None.
#    
# hclustplot1.old implements the case where couleur is a vector, hclustplotn
# the case of a matrix; hclustplot1 is a common wrapper.

if (exists("hclustplot1")) rm(hclustplot1);
hclustplot1 = function(hier1, couleur, title1="Colors sorted by hierarchical clustering", 
                     RowLabels=NULL, cex.RowLabels = 0.9, ...)
{
  if (length(dim(couleur))>1)
  {
    hclustplotn(hier1, couleur, RowLabels, cex.RowLabels, main = title1, ...);
  } else
  { 
    hclustplot1.old(hier1,couleur,title1);
  }
}

if (exists("hclustplot1.old")) rm(hclustplot1.old);
hclustplot1.old=function(hier1,couleur,title1="Colors sorted by hierarchical clustering") {
if (length(hier1$order) != length(couleur) ) {print("ERROR: length of color vector not compatible with no. of objects in the hierarchical tree")};
if (length(hier1$order) == length(couleur) ) {
barplot(height=rep(1, length(couleur)), col= as.character(couleur[hier1$order]),border=F, main=title1,space=0, axes=F)}
}
 
 
if (exists("hclustplotn")) rm(hclustplotn);
hclustplotn=function(hier1, Color, RowLabels=NULL, cex.RowLabels = 0.9, ...) 
{
  options(stringsAsFactors=FALSE);
  if (length(hier1$order) != dim(Color)[1] ) 
  { 
    stop("ERROR: length of color vector not compatible with no. of objects in the hierarchical tree");
  } else {
       No.Sets = dim(Color)[[2]];
       C = Color[hier1$order, ]; 
       step = 1/dim(Color)[[1]];
       ystep = 1/No.Sets;
       barplot(height=1, col = "white", border=F,space=0, axes=F, ...)
       for (j in 1:No.Sets)
       {
         for (i in 1:dim(C)[[1]])
         { 
           #lines(x=rep((i*step), times=2), y=c(ystep*(j-1),ystep*j),  col = as.character(C[i,j]));
           rect(xleft=((i-1)*step), ybottom=ystep*(j-1), xright = (i) * step, ytop = ystep*j,  
                    col = as.character(C[i,j]), border = as.character(C[i,j]));
         } 
         if (is.null(RowLabels))
         {
             text(as.character(j), pos=2, x=0, y=ystep*(j-0.5), cex=cex.RowLabels, xpd = TRUE);
         } else {
             text(RowLabels[j], pos=2, x=0, y=ystep*(j-0.5), cex=cex.RowLabels, xpd = TRUE);
         }
       }
       for (j in 1:No.Sets) lines(x=c(0,1), y=c(ystep*j,ystep*j));
  }
}
 




# ===================================================
# The function TOMplot1 creates a TOM plot
# Inputs:  distance measure, hierarchical (hclust) object, color label=couleur
if (exists("TOMplot1")) rm(TOMplot1);
TOMplot1=function(disttom,hier1, couleur,terrainColors=FALSE) {
no.nodes=length(couleur)
if (no.nodes != dim(disttom)[[1]] ) {print("ERROR: number of color labels does not equal number of nodes in disttom")} else {
labeltree=as.character(couleur)
labelrow  = labeltree
labelrow[hier1$order[length(labeltree):1]]=labelrow[hier1$order]
options(expressions = 10000)
if (terrainColors) heatmap(as.matrix(disttom),Rowv=as.dendrogram(hier1),Colv= as.dendrogram(hier1), scale="none",revC=T, ColSideColors=as.character(labeltree),RowSideColors=as.character(labelrow),
labRow=F, labCol=F, col = terrain.colors(1000)) else heatmap(as.matrix(disttom),Rowv=as.dendrogram(hier1),Colv= as.dendrogram(hier1), scale="none",revC=T, ColSideColors=as.character(labeltree),RowSideColors=as.character(labelrow),
labRow=F, labCol=F)
}
} #end of function


# ===================================================
# The function TOMplot2 creates a TOM plot where the top and left color bars can be different
# Inputs:  distance measure, hierarchical (hclust) object, color label=couleurTop, couleurLeft
if (exists("TOMplot2")) rm(TOMplot2);
TOMplot2=function(disttom,hier1, couleurTop, couleurLeft) {
no.nodes=length(couleurTop)
if (no.nodes != length(couleurLeft)) {stop("ERROR: number of top color labels does not equal number of left color labels")}
if (no.nodes != dim(disttom)[[1]] ) {stop("ERROR: number of color labels does not equal number of nodes in disttom")} else {
labeltree = as.character(couleurTop)
labelrow  = as.character(couleurLeft)
labelrow[hier1$order[length(labeltree):1]]=labelrow[hier1$order]
options(expressions = 10000)
heatmap(as.matrix(disttom),Rowv=as.dendrogram(hier1),Colv= as.dendrogram(hier1), scale="none", revC=T, ColSideColors=as.character(labeltree),RowSideColors=as.character(labelrow),
labRow=F, labCol=F)
}
} #end of function



# IGNORE THIS FUNCTION...
# The function "BestHeightCut" allows one to find the best height cut-off
# for a hierarchical clustering tree when external gene information is available
# It computes a Kruskal Wallis-test p-value for each height cut-off
# based on determining whether gene significance differs across branch membership.
if(exists("BestHeightCut")) rm(BestHeightCut);
BestHeightCut=function(hier1, GeneSignif, hcut=seq(0.1,.95,by=0.01) ) {
pvalues=rep(NA, length(hcut))
for (i in c(1:length(hcut))) {
colorhelp=modulecolor2(hier1,hcut[i])
if (length(unique(colorhelp))>1 ) {pvalues[i]=kruskal.test(GeneSignif, colorhelp)$p.value}
data.frame(hcut,pvalues)
}}




#####################################################################################################
################################################################################################################################
# D) Summing up modules using their first principal components (first eigengene)
#####################################################################################################
################################################################################################################################

# ===================================================
#The function ModulePrinComps1 finds the first principal component (eigengene) in each 
# module defined by the colors of the input vector "couleur" (Pardon my French).
# It also reports the variances explained by the first 5 principal components.
# This requires the R library impute
# The output is a list with 2 components: 
# 1) a data frame of module eigengenes (PCs), 
# 2) a data frame that lists the percent variance explained by the first 5 PCs of a module
# Options: if removeGrey=T, then no output is generated for the grey module.
# Recall that grey often denotes genes outside proper modules. 
if (exists("ModulePrinComps1" ) ) rm(ModulePrinComps1);
ModulePrinComps1=function(datexpr,couleur,removeGrey=F, FiveComponents=F) {
modlevels=levels(factor(couleur))
if ( removeGrey ) modlevels=setdiff(modlevels, c("grey") ); 
if (FiveComponents ) {print("To speed up the calculation, we only compute the five principal components of each module.  Therefore, the estimate of the proportion of variance explained is no longer accurate. If you want an accurate estimate of the proportion of var explained, please choose  the option FiveComponents=F ")  ;} 
PrinComps=data.frame(matrix(666,nrow=dim(datexpr)[[1]],ncol= length(modlevels))) 
varexplained= data.frame(matrix(666,nrow= 5,ncol= length(modlevels)))
names(PrinComps)=paste("PC",modlevels,sep="")
for(i in c(1:length(modlevels)) ){
print(i)   
modulename    = modlevels[i]
 restrict1= as.character(couleur)== modulename
 # in the following, rows are genes and columns are samples
    datModule=t(datexpr[, restrict1])
    saved.seed = .Random.seed;
    datModule=impute.knn(as.matrix(datModule))
    .Random.seed = saved.seed;
datModule=t(scale(t(datModule)))
if (FiveComponents ) { svd1 =svd(datModule, nu = 5, nv = 5) } else {svd1=svd(datModule)}
    mtitle=paste("PCs of ", modulename," module", sep="")
varexplained[,i]= (svd1$d[1:5])^2/sum(svd1$d^2)
# this is the first principal component
    pc1=svd1$v[,1]
 signh1=sign(sum(cor(pc1,  t(datModule))))
 if (signh1 != 0)  pc1=signh1* pc1
PrinComps[,i]= pc1
}
list(PrinComps=PrinComps, varexplained=varexplained)
}





#####################################################################################################
################################################################################################################################
# E) Relating a measure of gene significance to the modules 
#####################################################################################################
################################################################################################################################

# ===================================================
# The function ModuleEnrichment1 creates a bar plot that shows whether modules are enriched with
# significant genes.
# More specifically, it reports the mean gene significance for each module.
# The gene significance can be a binary variable or a quantitative variable. 
# It also plots the 95% confidence interval of the mean (CI=mean +/- 1.96* standard error).
# It also reports a Kruskal Wallis P-value.
if( exists("ModuleEnrichment1") ) rm(ModuleEnrichment1);
ModuleEnrichment1=function(genesignif1,couleur,title1="gene significance across modules",labely="Gene Significance",boxplot=F) {
if (length(genesignif1) != length(couleur) ) print("Error: vectors don't have the same lengths") else {
if (boxplot != TRUE) {
mean1=function(x) mean(x,na.rm=T) 
means1=as.vector(tapply(genesignif1,couleur,mean1));
se1= as.vector(tapply(genesignif1,couleur,stderr1))
par(mfrow=c(1,1))
barplot(means1,
names.arg=names(table(couleur) ),col= names(table(couleur) )
,ylab=labely)
err.bp(as.vector(means1), as.vector(1.96*se1), two.side=T)} else {
boxplot(split(genesignif1,couleur),notch=T,varwidth=T, col= names(table(couleur) ),ylab=labely)}

title(paste(title1,", p-value=", signif(kruskal.test(genesignif1,factor(couleur))$p.value,2)))
}
} # end of function


# IGNORE THIS...
# ===================================================
#The function fisherPvector allows one to compute Fisher exact p-values
# Thus it allows one to carry out an EASE analysis
# Output: a table of Fisher¡¯s exact p-values
# Input: annotation1 is a vector of gene annotations
# Input: couleur (French for color) denotes the module color of each gene
# Only those gene functions (levels of annotation1) that occur a certain mininum number of times
# (parameter= minNumberAnnotation) in the data set will be considered.  
if (exists("fisherPvector" ) ) rm(fisherPvector);
fisherPvector=function(couleur,annotation1,minNumberAnnotation=50) {
levelsannotation1=levels(annotation1)
levelscouleur=levels(factor(couleur))
no.couleur=length(levelscouleur)
restannotation1=table(annotation1)>minNumberAnnotation
no.annotation=sum( restannotation1)
datoutP=data.frame(matrix(666,nrow=no.annotation,ncol=no.couleur) )
#datoutProp=data.frame(matrix(666,nrow=no.annotation,ncol=2*no.couleur) )
#names(datoutProp)=paste("Prop",paste( rep(levelscouleur ,rep(2, length(levelscouleur))) ) , c("Y","N")  ,sep=".")
datoutProp=data.frame(matrix(666,nrow=no.annotation,ncol=no.couleur) )
names(datoutProp)=paste("Perc",levelscouleur , sep=".")
names(datoutP)=paste("p",levelscouleur,sep=".")
restlevelsannotation1= levelsannotation1[restannotation1]
row.names(datoutP)= restlevelsannotation1
for (i in c(1:no.annotation) ) {
for (j in c(1:no.couleur) ){
tab1=table( annotation1 !=restlevelsannotation1[i], couleur !=levelscouleur[j])
datoutP[i,j]=signif(fisher.test(tab1)$p.value,2) 
#datoutProp[i,2*j-1]=signif(tab1[1,1]/sum(tab1[,1] ),2)
#datoutProp[i,2*j]= signif(tab1[1,2]/sum(tab1[,2]) ,2)
} 
table2=table(annotation1 !=restlevelsannotation1[i], couleur)
datoutProp[i,]= signif(table2[1,]/apply(table2,2,sum),2)
}
data.frame(datoutP,datoutProp)
} # end of function fisherPvector



#####################################################################################################
################################################################################################################################
# F) Carrying out a within module analysis (computing intramodular connectivity etc) 
#####################################################################################################
################################################################################################################################

# ===================================================
#The function DegreeInOut computes for each gene 
#a) the total number of connections, 
#b) the number of connections with genes within its module, 
#c) the number of connections with genes outside its module
# When scaledToOne=TRUE, the within module connectivities are scaled to 1, i.e. the max(K.Within)=1 for each module
if (exists("DegreeInOut")) rm(DegreeInOut); DegreeInOut =function(adj1, couleur,scaledToOne=FALSE) {
no.nodes=length(couleur)
couleurlevels=levels(factor(couleur))
no.levels=length(couleurlevels)
kWithin=rep(-666,no.nodes )
diag(adj1)=0
for (i in c(1:no.levels) ) {
rest1=couleur==couleurlevels[i];
if (sum(rest1) <3 ) { kWithin[rest1]=0 } else {
kWithin[rest1]=apply(adj1[rest1,rest1],2,sum)
if (scaledToOne) kWithin[rest1]=kWithin[rest1]/max(kWithin[rest1])}
}
kTotal= apply(adj1,2,sum) 
kOut=kTotal-kWithin
if (scaledToOne) kOut=NA
kDiff=kWithin-kOut
data.frame(kTotal,kWithin,kOut,kDiff)
}


# =======================================================================
# The function WithinModuleCindex1 relates the node measures (e.g. connectivities) 
# to "external" node significance information within each  module,
# i.e. it  carries out a by module analysis.
# Output: first column reports the spearman correlation p-value between the network variable and the 
# node significance. The next columns contain the Spearman correlations between the variables.
if (exists("WithinModuleAnalysis1")) rm(WithinModuleAnalysis1);
WithinModuleAnalysis1=function(datnetwork,nodesignif, couleur) {
cortesthelp=function( x ) {
len1=dim(x)[[2]]-1
out1=rep(666, len1);
for (i in c(1:len1) ) {out1[i]= signif( cor.test(x[,i+1], x[,1], method="s",use="p" )$p.value ,2) }
data.frame( variable=names(x)[-1] , NS.CorPval=out1, NS.cor=t(signif(cor (x[,1], x[,-1],use="p",method="s"),2)), signif(cor(x[,-1],use="p",method="s"),2) )
} #end of function cortesthelp
print("IGNORE  the warnings...");
by( data.frame(nodesignif, datnetwork),couleur,cortesthelp);
} #end of function WithinModuleAnalysis




# =======================================================================
# The function WithinModuleCindex1 relates the node measures (e.g. connectivities) 
# to "external" node significance information within each  module, 
# i.e. it  carries out a by module analysis.
# BUT it focuses on the C-index also known as area under the ROC curve
# This measure is related to Kendall's Tau statistic and Somer's D, 
# see F. Harrel (Regression Modeling Strategies). Springer. 
# It requires the following library
library(Hmisc)
# Output: the first column reports the C-index and the second, p-value 
if (exists("WithinModuleCindex1")) rm(WithinModuleCindex1);
WithinModuleCindex1=function(datnetwork,nodesignif, couleur) {
CindexFun=function( x ) {
len1=dim(x)[[2]]-1
outC=rep(666, len1);
outP=rep(666, len1);
for (i in c(1:len1) ) {rcor1=rcorr.cens(x[,i+1], x[,1])
outC[i]=rcor1[[1]] 
outP[i]=1- pnorm(abs(rcor1[[2]]/rcor1[[3]]))
}
data.frame( variable=names(x)[-1] , C.index=outC, p.value=outP)
} #end of function CindexFun
#print("IGNORE  the warnings...");
by( data.frame(nodesignif, datnetwork),couleur,CindexFun);
} #end of function WithinModuleAnalysis


# The following function allows on to plot a gene (node) significance measure versus
# connectivity.
if(exists("plotConnectivityGeneSignif1") ) rm( plotConnectivityGeneSignif1);
plotConnectivityGeneSignif1=function(degree1,genesignif1,color1="black", 
title1="Gene Significance vs Connectivity" , xlab1="Connectivity", ylab1="GeneSignificance") {
lm1=lm(genesignif1~degree1 ,na.action="na.omit")
plot(degree1, genesignif1, col=color1,ylab=ylab1,xlab=xlab1,main=paste(title1, ", cor=",  
signif(cor( genesignif1,degree1, method="s",use="p" )   ,2) ))
abline(lm1)
}




 var1=function(x) var(x, na.rm=T)
 no.present=function(x) sum(!is.na(x))
 




# The following function computes the module eigengene based connectivity.
# Input: datExpr= a possibly very large gene expression data set (columns represent genes)
# datPC=data frame of module eigengenes (columns correspond to module eigengenes or PCs)
# A module eigengene based connectivity KME value will be computed if the gene has 
# a non missing expression value in at least MinimumNoSamples arrays.
# Output a data frame where columns are the KME values corresponding to different modules.
# By splitting the expression data into different batches, the function can deal with 
# tens of thousands of gene expression data. 
# If there are many eigengenes (say hundreds) consider decreasing the batch size.
if (exists("signedKME") ) rm(signedKME); 
signedKME=function(datExpr, datPC, batchsize=20000, MinimumNoSamples=5) {
 datExpr=data.frame(datExpr)
 datPC=data.frame(datPC)
 datKME=list()
 if (dim(datPC)[[1]] != dim(datExpr)[[1]] ) print("ERROR in signedKME function: dim(datPC)[[1]] != dim(datExpr)[[1]] ")
 varianceZeroIndicatordatExpr=as.vector(apply(datExpr,2,var1))==0
 varianceZeroIndicatordatPC =as.vector(apply(datPC,2,var1))==0
 if (sum(varianceZeroIndicatordatExpr,na.rm=T)>0 ) print("Warning: in signedKME: remove constant columns (zero variance) in datExpr" )
 if (sum(varianceZeroIndicatordatPC,na.rm=T)>0 ) print("ERROR in signedKME: remove constant columns in datPC" )
 no.presentdatExpr=as.vector(apply(datExpr,2, no.present) )
 if (min(no.presentdatExpr)<5 ) print("Warning for signedKME: some gene expressions have fewer than 5 observations. Please remove genes with too many missing values")
 if (dim(datPC)[[1]] == dim(datExpr)[[1]] & sum(varianceZeroIndicatordatPC,na.rm=T)==0 ) {
no.genes=dim(datExpr)[[2]]
no.samples=dim(datExpr)[[1]]
no.batches=as.integer(no.genes/ batchsize)
 datKME=data.frame(matrix(NA, nrow=dim(datExpr)[[2]], ncol=dim(datPC)[[2]]))
 names(datKME)=paste("kME", substring(names(datPC), first=3, last=100), sep="")  
 dimnames(datKME)[[1]] = names(datExpr) 
if (no.batches>0) {
for (i in 1:no.batches) {
print(paste("batch number = ", i))
index1=c(1:batchsize)+(i-1)* batchsize
datKME[index1,]= cor(datExpr[,index1], datPC ,use="p")     
# If fewer than MinimumNoSamples contain gene expression information for a given
# gene, then we set its connectivity to 0.
NoSamplesAvailable=apply(!is.na(datExpr[,index1]),2,sum)
datKME[index1,][NoSamplesAvailable< MinimumNoSamples,]=NA
} # end of for (i in 1:no.batches
} # end of if (no.batches>0)...
if (no.genes-no.batches*batchsize>0 ) {
restindex=c((no.batches*batchsize+1):no.genes)
datKME[restindex,] = cor(datExpr[,restindex], datPC ,use="p")     
NoSamplesAvailable=apply(!is.na(datExpr[,restindex]),2,sum)
datKME[restindex,][NoSamplesAvailable< MinimumNoSamples,]=NA
} # end of if
 } # end of if statement checking input data
 datKME
 } # end of function signedKME





#####################################################################################################
################################################################################################################################
# G) Miscellaneous other functions, e.g. for computing the cluster coefficient.
#####################################################################################################
################################################################################################################################



# ===================================================
# The function ClusterCoef.fun computes the cluster coefficients.
# Input is an adjacency matrix 
if(exists("ClusterCoef.fun")) rm(ClusterCoef.fun) ; ClusterCoef.fun=function(adjmat1) {
diag(adjmat1)=0
no.nodes=dim(adjmat1)[[1]]
computeLinksInNeighbors <- function(x, imatrix){x %*% imatrix %*% x}
nolinksNeighbors <- c(rep(-666,no.nodes))
total.edge <- c(rep(-666,no.nodes))
maxh1=max(as.dist(adjmat1) ); minh1=min(as.dist(adjmat1) ); 
if (maxh1>1 | minh1 < 0 ) {print(paste("ERROR: the adjacency matrix contains entries that are larger than 1 or smaller than 0!!!, max=",maxh1,", min=",minh1)) } else { 
nolinksNeighbors <- apply(adjmat1, 1, computeLinksInNeighbors, imatrix=adjmat1)
plainsum  <- apply(adjmat1, 1, sum)
squaresum <- apply(adjmat1^2, 1, sum)
total.edge = plainsum^2 - squaresum
CChelp=rep(-666, no.nodes)
CChelp=ifelse(total.edge==0,0, nolinksNeighbors/total.edge)
CChelp}
} # end of function



# ===================================================
# The function err.bp  is used to create error bars in a barplot
# usage: err.bp(as.vector(means), as.vector(stderrs), two.side=F)
err.bp<-function(daten,error,two.side=F){
 if(!is.numeric(daten)) {
      stop("All arguments must be numeric")}
 if(is.vector(daten)){ 
    xval<-(cumsum(c(0.7,rep(1.2,length(daten)-1)))) 
 }else{
    if (is.matrix(daten)){
      xval<-cumsum(array(c(1,rep(0,dim(daten)[1]-1)),
dim=c(1,length(daten))))+0:(length(daten)-1)+.5
    }else{
      stop("First argument must either be a vector or a matrix") }
 }
 MW<-0.25*(max(xval)/length(xval)) 
 ERR1<-daten+error 
 ERR2<-daten-error
 for(i in 1:length(daten)){
    segments(xval[i],daten[i],xval[i],ERR1[i])
    segments(xval[i]-MW,ERR1[i],xval[i]+MW,ERR1[i])
    if(two.side){
      segments(xval[i],daten[i],xval[i],ERR2[i])
      segments(xval[i]-MW,ERR2[i],xval[i]+MW,ERR2[i])
    } 
 } 
} 

# ===================================================
# this function computes the standard error
if (exists("stderr1")) rm(stderr1)
stderr1 <- function(x){ sqrt( var(x,na.rm=T)/sum(!is.na(x))   ) }



# ===================================================
# The following two functions are for displaying the pair-wise correlation in a panel when using the command "pairs()"
# Typically, we use "pairs(DATA, upper.panel=panel.smooth, lower.panel=panel.cor, diag.panel=panel.hist)" to
# put the correlation coefficients on the lower panel.
panel.cor <- function(x, y, digits=2, prefix="", cex.cor){
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex * r)
}
panel.hist <- function(x, ...){
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}



# ===================================================
# this function computes the standard error
if (exists("stderr1")) rm(stderr1);
stderr1 <- function(x){ sqrt( var(x,na.rm=T)/sum(!is.na(x))   ) }




# ===================================================
# This function collects garbage
if (exists("collect_garbage")) rm(collect_garbage);
collect_garbage=function(){while (gc()[2,4] != gc()[2,4] | gc()[1,4] != gc()[1,4]){}}
collect_garbage()


# this function is used for computing the Rand index below...
# ===================================================
if (exists("choosenew") ) rm(choosenew)
choosenew <- function(n,k){
  n <- c(n)
  out1 <- rep(0,length(n))
  for (i in c(1:length(n)) ){
    if (n[i]<k) {out1[i] <- 0}
    else {out1[i] <- choose(n[i], k)}}
  out1	
}


# ===================================================
# the following function computes the Rand index between 2 clusterings
# assesses how similar two clusterings are
if (exists("Rand1") ) rm(Rand1)
Rand2 <- function(tab,adjust=T) {
  a <- 0; b <- 0; c <- 0; d <- 0; nn <- 0
  m <- nrow(tab);
  n <- ncol(tab);
  for (i in 1:m) {
    c<-0
    for(j in 1:n) {
      a <- a+choosenew(tab[i,j],2)
      nj <- sum(tab[,j])
      c <- c+choosenew(nj,2)
    }
    ni <- sum(tab[i,])
    b <- b+choosenew(ni,2)
    nn <- nn+ni
  }
  if(adjust==T) {
    d <- choosenew(nn,2)
    adrand <- (a-(b*c)/d)/(0.5*(b+c)-(b*c)/d)
    adrand
  } else {
    b <- b-a
    c <- c-a
    d <- choosenew(nn,2)-a-b-c
    rand <- (a+d)/(a+b+c+d)
    rand
  }
}

# ===================================================
# This function is used in "pairs()" function. The problem of the original  panel.cor is that 
# when the correlation coefficient is very small, the lower panel will have a large font 
# instead of a mini-font in a saved .ps file. This new function uses a format for corr=0.2 
# when corr<0.2, but it still reports the original value of corr, with a minimum format.

panel.cor1=function(x, y, digits=2, prefix="", cex.cor){
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    txt1=txt
    r1=r
    if (r<0.2) {
        r1=0.2
        txt1 <- format(c(r1, 0.123456789), digits=digits)[1]
        txt1 <- paste(prefix, txt1, sep="")
        }
    if(missing(cex.cor)) cex <- 0.8/strwidth(txt1)
    cex = cex * r1
    r <- round(r, digits)
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    text(0.5, 0.5, txt, cex=cex)
}

 







#merge the minor cluster into the major cluster
merge2Clusters = function(mycolorcode, mainclusterColor, minorclusterColor){
  mycolorcode2 = ifelse(as.character(mycolorcode)==minorclusterColor, mainclusterColor, as.character(mycolorcode) )
  fcolorcode   =factor(mycolorcode2)
  fcolorcode
}



# A set of global variables and functions that should help handling color names for some 400+ modules.
# A vector called GlobalStandardColors is defined that holds color names with first few entries 
# being the well-known and -loved colors. The rest is randomly chosen from the color names of R,
# excluding grey colors.
#---------------------------------------------------------------------------------------------------------
## GlobalStandardColors 
#
#---------------------------------------------------------------------------------------------------------
# This code forms a vector of color names in which the first entries are given by BaseColors and the rest
# is "randomly" chosen from the rest of R color names that do not contain "grey" nor "gray".

BaseColors = c("turquoise","blue","brown","yellow","green","red","black","pink","magenta",
"purple","greenyellow","tan","salmon","cyan", "midnightblue", "lightcyan","grey60", "lightgreen", 
"lightyellow", "royalblue", "darkred", "darkgreen", "darkturquoise", "darkgrey", "orange", "darkorange",
 "white", "skyblue", "saddlebrown", "steelblue", "paleturquoise", "violet", "darkolivegreen",
 "darkmagenta", "white" );
RColors = colors()[-grep("grey", colors())];
RColors = RColors[-grep("gray", RColors)];
InBase = match(BaseColors, RColors);
ExtraColors = RColors[-c(InBase[!is.na(InBase)])];
No.Extras = length(ExtraColors);
# Here is the vector of colors that should be used by all functions:
GlobalStandardColors = c(BaseColors, ExtraColors[rank(sin(13*c(1:No.Extras) +sin(13*c(1:No.Extras))) )] );
rm(BaseColors, RColors, ExtraColors, No.Extras);
#---------------------------------------------------------------------------------------------------------
#
# ColorsFromLabels
#
#---------------------------------------------------------------------------------------------------------
# This function converts integer numerical labels Labels into color names in the order either given by
# StandardColors,
# or (if StandardColors==NULL) by GlobalStandardColors. If GreyIsZero == TRUE, labels 0 will be assigned
# the color grey; otherwise presence of labels below 1 will trigger an error.
ColorsFromLabels = function(Labels, ZeroIsGrey = TRUE, StandardColors = NULL)
{
  if (is.null(StandardColors)) StandardColors = GlobalStandardColors;
  if (ZeroIsGrey) MinLabel = 0 else MinLabel = 1
  if (sum( (Labels>=MinLabel) & (Labels <= length(StandardColors)) )!= length(Labels))
    stop(paste("Input error: something's wrong with Labels. Either they are not a numeric vector,",
               "or some values are below", MinLabel,  "or some values are above the maximum number of colors", length(StandardColors)));
  Colors = rep("grey", length(Labels));
  Colors[Labels!=0] = StandardColors[Labels[Labels!=0]];
  Colors;}
#---------------------------------------------------------------------------------------------------------
#
# DisplayColors
#
#---------------------------------------------------------------------------------------------------------
# This function plots a barplot with all colors given. If Colors are not given, GlobalStandardColors are
# used, i.e. if you want to see the GlobalStandardColors, just call this function without parameters.
DisplayColors = function(Colors = NULL)
{
  if (is.null(Colors)) Colors = GlobalStandardColors;
  barplot(rep(1, length(Colors)), col = Colors, border = Colors);}


#----------------------------------------------------------## ColorsFromLabels##------------------------------



###############################################################################
# I) Functions for merging modules based on a high correlation of the module eigengenes
###############################################################################
# The following network functions were created by Peter Langfelder.


# This function is an alternative to the print command that enforces immediate output.

if (exists("print.flush")) { remove(print.flush); collect_garbage(); }
print.flush = function(...)
{
  x = print(...)
  if (exists("flush.console")) x=flush.console();
}

if (exists("PrintSpaces")) { remove(PrintSpaces); collect_garbage(); }
PrintSpaces = function(print.level)
{
  if (print.level>0) 
    {
      spaces = paste(rep("  ", times=print.level), collapse="");
    } else
    {
      spaces = "";
    }
  spaces;
}




# This function forms the average gene expression for each color in the color vector.
# Input: gene expression data (rows are samples, columns are genes)
# In many applications, we recommend to input normalized gene expression values, i.e.
# NormExprData=scale(datExpr). colors= is a vector that encodes the module color of each gene.

AverageExprMatrix = function(NormExprData, colors) {
  no.genes = dim(NormExprData)[2]
  no.samples = dim(NormExprData)[1]
  colorsf = as.factor(colors)
  AverageExpr = matrix(ncol=nlevels(colorsf), nrow = no.samples)
  ExprDataMtrx = as.matrix(NormExprData)
  for (i in (1:no.samples)) AverageExpr[i,] = tapply(ExprDataMtrx[i,], colorsf, mean)
  AverageExpr
}

#---------------------------------------------------------------------------------
#
# ModuleNumber. Function from Peter Langfelder
#
#---------------------------------------------------------------------------------
# Similar to modulecolor2 above, but returns numbers instead of colors, which is oftentimes more useful.
# 0 means unassigned.
# This is particularly useful when dealing with a lot of branches in the tree.
# Or when integer values instead of strings are preferred to label modules.
# Return value is a simple vector, not a factor.
# Caution: the module numbers are neither sorted nor sequential, the only guarranteed fact is that grey
# probes are labeled by 0 and all probes belonging to the same module have the same number.


ModuleNumber = function(HierTree, CutHeight = 0.9, MinSize = 50)
{
  Branches = cutree(HierTree, h = CutHeight);
  NOnBranches = table(Branches);
  #NOnBranches[i][[1]] somehow gives the number of elements on branch i.
  TrueBranch = NOnBranches >= MinSize;
  Branches[!TrueBranch[Branches]] = 0;
  
  #NewLabels = levels(factor(Branches));
  #for (lab in 1:length(NewLabels)) if (NewLabels[lab]!=0)
  #  Branches[Branches==NewLabels[lab]] = lab;

  Branches;

}

#-------------------------------------------------------------------------------------
#
#  ModulePrincipalComponents
#
#-------------------------------------------------------------------------------------
# This function is an alternative to ModulePrincomps1. Probably it produces the same output.
# This function computes the principal components of modules of a given network.
# The first PC of each module is the module eigengene.
# Input: Data: expression data, module colors. 
# In general, the first PC is only defined up to a sign.
# This function optionally orients the first PC so that it has a positive correlation
# with the average normalized gene expression profile, see the function AverageExpression
# AlignPCs can take the values "", "along average".
# output : a dataframe of module eigengenes (PCs)

ModulePrincipalComponents = function(Data, ModuleColors, AlignPCs = "along average", verbose = 1,
  print.level=0) 
{
  spaces = PrintSpaces(print.level);

  AlignPCsRecognizedValues =  c("", "along average");
  if (!is.element(AlignPCs, AlignPCsRecognizedValues)) {
    print.flush(paste("ModulePrincipalComponents: Error:",
                "parameter AlignPCs has an unrecognised value:", 
                AlignPCs, "; Recognized values are ", AlignPCsRecognizedValues));
    stop()
  }

  if (verbose>0) print.flush(paste(spaces, "ModulePrincipalComponents: Calculating PCs"));
  FullPCs = ModulePrinComps1(Data, ModuleColors);
  PCs = FullPCs$PrinComps;

  if (AlignPCs == "") AlignedPCs = PCs;
  if (AlignPCs == "along average") 
  {
    if (verbose>0) print.flush(paste(spaces,"ModulePrincipalComponents:", 
                     "Aligning PCs with average expression for each module."))
    if (verbose>1) print.flush(paste(spaces,"  ++ Calculating averages..."));
    NormData = scale(Data);
    AverageModuleExpr = data.frame(AverageExprMatrix(NormData, ModuleColors));
    if (verbose>1) print.flush(paste(spaces,"  ++ Aligning principal components..."));
    AverageAndPCCor = diag(cor(PCs, AverageModuleExpr, use = "p"));
    sign.matrix = matrix(0, nrow = dim(PCs)[2], ncol = dim(PCs)[2]);
    diag(sign.matrix) = sign(AverageAndPCCor);
    AlignedPCs = as.data.frame(as.matrix(PCs) %*% sign.matrix);
    names(AlignedPCs) = names(PCs);
    rownames(AlignedPCs) = rownames(PCs);
    names(AverageModuleExpr) = names(PCs);
    rownames(AverageModuleExpr) = rownames(PCs);
    if (verbose>1) print.flush(paste(spaces,"  ++ done."));
  }
  RetPCs = list(data = AlignedPCs, VarExplained = FullPCs$varexplained, 
                ModuleConformity = FullPCs$ModuleConformity, AverageExpr = AverageModuleExpr);
  RetPCs;
}


#--------------------------------------------------------------------------------------
# 
# ModulePCs
# This helper function is used for the function MergeCloseModules.
#--------------------------------------------------------------------------------------
# Input: Data are gene expression data where columns correspond to genes 
# and rows correspond to samples that may correspond to different subsets (e.g. study 1 and study 2)
# Subsets is a vector of labels "1", "2", "3", etc (i.e. positive integers).
# Output: a vector of lists where each list contains the following
# data=module eigengenes, 
# AverageExpr = average normalized module expression
# VarExplained= variance explained by the eigengenes 
# Options: AlignPCs= chooses the sign of each module eigengene by
# enforcing a positive correlation with the average normalized expression.


ModulePCs = function(Data, Subsets, ModuleColors, OnlySet = NULL,
                     AlignPCs="along average", verbose=1, print.level=0)
{
  spaces = PrintSpaces(print.level)
  SubsetFactor = factor(Subsets);
  No.Subsets = nlevels(SubsetFactor);
  if (verbose>0) print.flush(paste(spaces,"ModulePCs: Looking for module PCs."));
  PCs = vector(mode="list", length=No.Subsets);
  if (is.null(OnlySet))
  {
      CalculatedSubsets = c(1:No.Subsets);
  } else {
      CalculatedSubsets = c(OnlySet);
  }
  for (subs.ind in CalculatedSubsets) {
    subset = levels(SubsetFactor)[subs.ind];
    if (verbose>0) print.flush(paste(spaces,"  Working on subset", as.character(subset), "...")); 
    SubsetNetworkData = Data[SubsetFactor == subset, ];
    SubsetColors = ModuleColors; 
    SubsetPCs = ModulePrincipalComponents(Data = SubsetNetworkData, 
                          ModuleColors = SubsetColors, AlignPCs = AlignPCs, verbose = verbose-1,
                          print.level = print.level+1);
    PCs[[subs.ind]] = list(data = SubsetPCs$data, AverageExpr = SubsetPCs$AverageExpr, 
                           ModuleConformity = SubsetPCs$ModuleConformity, 
                           VarExplained = SubsetPCs$VarExplained);
    rm(SubsetNetworkData); rm(SubsetColors); rm(SubsetPCs); collect_garbage();
  }
  PCs;
}


#---------------------------------------------------------------------------------------------
#
# MergeCloseModules
#
#---------------------------------------------------------------------------------------------
# This function merges modules whose PCs fall on one branch of a hierarchical clustering tree
# Method: 
# First, the module eigenes are computed corresponding to each color (but see option PCs)
# Second, we use average linkage hierarchical clustering of the module eigengenes to arrive at a dendrogram
# Third, branches are cut-off the dendrogram using a given height (CutHeight)
# Fourth, modules whose PCs fall on one branch are merged.
# If the option Relabel=True, then the colors of the resulting merged modules 
# are assigned according to size, e.g. turquoise encodes the largest module
# If Relabel=F, then the color of the merged module is chosen according to that of the first module.
# Options: the parameter PCs can be used to input a vector of lists 
# where each list must have a component named data that contains the data frame of module eigengenes.
# UseAbs specifies whether the absolute value of the correlation should be used for constructing the clustering tree of the PCs.
# if CalculateNewPCs = TRUE is specified then the PCs are newly calculated based on the new merged colors.
# The options Subsets=NULL, OnlySet = NULL are useful when the data are comprised of multiple subsets (e.g. study 1 and 2).

# StandardColors allows one to specify your own order of colors to be assigned according to module size.
# Thus, specifying StandardColors only makes sense if Relabel=True.


if(exists("MergeCloseModules")) rm(MergeCloseModules);
MergeCloseModules = function(datE, OriginalColors, CutHeight=0.20, 
                             PCs = NULL, UseAbs = F, 
                             Relabel = FALSE, StandardColors = NULL, CalculateNewPCs = TRUE,
                             Subsets = NULL, OnlySet = NULL, 
                             verbose = 1, print.level=0)
{

  PCsInSingleFrame = FALSE;
  spaces = PrintSpaces(print.level);

  if (verbose>0) print.flush(paste(spaces, 
            "MergeCloseModules: Merging modules whose distance is less than", CutHeight));

  if (is.null(Subsets)) 
  {
    Subsets = rep(1, dim(datE)[1]);
    if (!is.null(PCs))
    {
      if (is.null(PCs[[1]]$data))
      {
        if (verbose>1) print.flush(paste(spaces, "  PCs appear to be a single data frame - will work on",
                                         "that assumption."));
        if (length(dim(PCs))!=2)
          stop(paste("PCs must be given either as a vector of lists, one list for each subset, with",
                     "element 'data' containing the frame or matrix of eigengenes",
                     " OR for a single data set PCs can be given as a dataframe or matrix of correct",
                     "dimensions."));
        PCsX = vector(mode="list", length = 1);
        PCsX[[1]] = list(data = PCs);
        rm(PCs); PCs = PCsX; rm(PCsX);
        PCsInSingleFrame = TRUE;
      } 
      if (dim(PCs[[1]]$data)[1]!=dim(datE)[1])
        stop(paste("Number of samples in PCs is incompatible with number of samples in datE."));
    }
  } else
  {
    SubsetsX = Subsets;
    Subsets = as.numeric(factor(Subsets));
  }

  if (!is.null(PCs))
  {
    SubsetsX = Subsets;
    No.Sets = max(Subsets);
    if (length(PCs)!=No.Sets)
      stop(paste("Number of sets given by Subsets is incompatible with the number of eigengene sets",
                 "given in PCs."));
    SamplesInSet = table(Subsets);
    for (i in 1:No.Sets)
    {
      if (dim(PCs[[i]]$data)[1]!=SamplesInSet[i])
          stop(paste("Number of samples in PCs is incompatible with subset length for subset",
                     levels(factor(SubsetsX))[i], " (first occurence)."));
    }
  }

  if (dim(datE)[1]!=length(Subsets))
    stop("Number of genes in datE is different from the length of given subset vector. They must equal.");

  if (dim(datE)[2]!=length(OriginalColors))
    stop("Number of genes in datE is different from the length of original colors. They must equal.");

  if ((CutHeight <0) | (CutHeight>(1+as.integer(UseAbs)))) 
    stop(paste("Given CutHeight is out of sensible range between 0 and", 1+as.integer(UseAbs) ));

  # If ordered PCs were not given, calculate them

  if (is.null(PCs)) 
  {
    PCs = ModulePCs(datE, Subsets, ModuleColors = OriginalColors,
                    OnlySet = OnlySet, 
                    verbose = verbose-1, print.level = print.level+1);
    collect_garbage();
  } else if (nlevels(as.factor(OriginalColors))!=dim(PCs[[1]]$data)[2])
  {
    if (verbose>0) print.flush(paste(spaces, "  Number of given module colors", 
              "does not match number of given MEs => recalculating the MEs."))
    PCs = ModulePCs(datE, Subsets, ModuleColors = OriginalColors,
                    OnlySet = OnlySet, 
                    verbose = verbose-1, print.level = print.level+1);
    collect_garbage();
  }

  # Cluster the found module eigengenes and merge ones that are too close to one another _in both sets_.

  No.Sets = nlevels(as.factor(Subsets));
  
  MEDiss = vector(mode="list", length = No.Sets);
  if (is.null(OnlySet))
  {
      CalculatedSubsets = c(1:No.Sets);
  } else {
      CalculatedSubsets = c(OnlySet);
  }
  for (set in CalculatedSubsets)
  {
    IndexRange = c(1:(nlevels(as.factor(OriginalColors))-1));
    if (UseAbs)
    {
        diss = 1-abs(cor(PCs[[set]]$data[, IndexRange], use = "p"));
    } else {
        diss = 1-cor(PCs[[set]]$data[, IndexRange], use = "p");
    }
    MEDiss[[set]] = list(Diss = diss);
  }
  
  if (is.null(OnlySet))
  {
    ConsDiss = (MEDiss[[1]]$Diss)
    if (No.Sets>1) for (set in 2:No.Sets)
      ConsDiss = pmax(ConsDiss, MEDiss[[set]]$Diss);
  } else {
    ConsDiss = MEDiss[[OnlySet]]$Diss;
  }
  
  METree = hclust(as.dist(ConsDiss), method = "average");
  METreeBranches = as.factor(ModuleNumber(HierTree = METree, CutHeight = CutHeight, MinSize = 1));

  # Analyze the branches: look for the ones that contain more than one original module
  
  MEUniqueBranches = levels(METreeBranches);
  MENo.Branches = nlevels(METreeBranches)
  MENumberOnBranch = rep(0, times = MENo.Branches);
  for (branch in 1:MENo.Branches)
  {
    MENumberOnBranch[branch] = sum(METreeBranches==MEUniqueBranches[branch]);
  }
  
  MergedColors = OriginalColors;
  
  # Merge modules on the same branch

  for (branch in 1:MENo.Branches) if (MENumberOnBranch[branch]>1)
  {
    ModulesOnThisBranch = names(METreeBranches)[METreeBranches==MEUniqueBranches[branch]];
    ColorsOnThisBranch = substring(ModulesOnThisBranch, 3);
    if (verbose>3) print.flush(paste(spaces, "  Merging original colors", paste(ColorsOnThisBranch, 
                                         collapse=", ")));
    for (color in 2:length(ColorsOnThisBranch))
      MergedColors[MergedColors==ColorsOnThisBranch[color]] = ColorsOnThisBranch[1];
  }
  
  No.Mods = nlevels(as.factor(MergedColors));
  RawModuleColors = levels(as.factor(MergedColors));

  # print(paste("No. of new modules: ", No.Mods));
  # print(paste("Merged module colors:"));
  # print(table(as.factor(MergedColors)));
 
  if (Relabel) 
  {
     # Relabel the merged colors to the usual order based on the number of genes in each module
 
     if (is.null(StandardColors))
     {
       StandardColors = c("turquoise","blue","brown","yellow","green","red","black","pink","magenta",
                   "purple","greenyellow","tan","salmon","cyan", "midnightblue", "lightcyan","grey60", 
                   "lightgreen", "lightyellow", "royalblue", "darkred", "darkgreen", "darkturquoise", 
                   "darkgrey", "orange", "darkorange", "white" );
   
     }
     
     No.GenesInModule = rep(0, No.Mods);
     for (mod in 1:No.Mods) No.GenesInModule[mod] = sum(MergedColors==RawModuleColors[mod]);
   
     # print(paste("No.GenesInModule: ", paste(No.GenesInModule, collapse=", ")));
     
     SortedRawModuleColors = RawModuleColors[order(-No.GenesInModule)]

     # print(paste("SortedRawModuleColors:", paste(SortedRawModuleColors, collapse=", ")));
     
     # Change the color names to the standard sequence, but leave grey grey (that's why rank in general does
     # not equal color)
 
     MergedNewColors = MergedColors;
 
     if (verbose>3) print(paste(spaces, "   Changing original colors:"));
     rank = 0;
     for (color in 1:length(SortedRawModuleColors)) if (SortedRawModuleColors[color]!="grey")
     {
       rank = rank + 1;
       if (verbose>3) print(paste(spaces, "      ", SortedRawModuleColors[color], 
                                  "to ", StandardColors[rank]));
       MergedNewColors[MergedColors==SortedRawModuleColors[color]] = StandardColors[rank];
     }
     # print("Table of new MergedColors:");
     # print(table(as.factor(MergedNewColors)));
  } else {
     MergedNewColors = MergedColors; 
  }

  if (CalculateNewPCs)
  {
     if (verbose>0) print.flush(paste(spaces, "  Calculating new PCs..."));
     NewPCs = ModulePCs(datE, Subsets, ModuleColors = MergedNewColors,
                        OnlySet = OnlySet, 
                        verbose = verbose-1, print.level = print.level+1);
  } else {
     NewPCs = NULL;
  }

  if (PCsInSingleFrame) 
  {
    NewPCs = NewPCs[[1]]$data;
    PCs = PCs[[1]]$data;
  }

  list(Colors = MergedNewColors, ClustTree = METree, CutHeight = CutHeight, OldPCs = PCs, NewPCs = NewPCs);
}





###############################################################################################################################
# I) GENERAL STATISTICAL FUNCTIONS

mean1=function(x) mean(x,na.rm=T)
var1=function(x) var(x,na.rm=T)



if(exists("scatterplot1") ) rm(scatterplot1);
scatterplot1=function(x,y, title1="",col1="black",xlab1=NA,ylab1=NA, cex1=1, cex.axis1=1.5,cex.lab1=1.5, cex.main1=1.5 ,ylim1=-1,correlationmethod="p" ){
if ( is.na(xlab1) ) xlab1= as.character(match.call(expand.dots = FALSE)$x)
if ( is.na(ylab1) ) ylab1= as.character(match.call(expand.dots = FALSE)$y)
x= as.numeric(as.character(x))
y= as.numeric(as.character(y))
cor1=signif(cor(x,y,use="p",method=correlationmethod),2)
corp=signif(cor.test(x,y,use="p",method=correlationmethod)$p.value,2)
if (corp<10^(-20) ) corp="<10^{-20}"
if ( length(ylim1)==2) {plot(x,y, main=paste(title1,"cor=", cor1,"p=",corp),col=col1,xlab=xlab1,ylab=ylab1, cex=cex1, 
cex.axis=cex.axis1,cex.lab=cex.lab1, cex.main=cex.main1,ylim=ylim1)} else {
plot(x,y, main=paste(title1,"cor=", cor1,"p=",corp),col=as.character(col1),xlab=xlab1,ylab=ylab1, cex=cex1, cex.axis=cex.axis1,cex.lab=cex.lab1, cex.main=cex.main1)
}
}



if  (exists("boxplot1") ) rm(boxplot1);
boxplot1=function(x,g,title1=" ",xlab1=NA,ylab1=NA, cex.axis1=1.5,cex.lab1=1.5, cex.main1=1.5,cex1=1.5,col1="white") {
if ( is.na(xlab1) ) xlab1= as.character(match.call(expand.dots = FALSE)$x)
print(xlab1)
if ( is.na(ylab1) ) ylab1= as.character( match.call(expand.dots = FALSE)$g)
print(ylab1)
p1=signif(kruskal.test(x, factor(g) )$p.value,2)
if (p1< 5.0*10^(-22) ) p1="<5.0x10^{-22}"
boxplot(x~factor(g), notch=T,varwidth=T, main=paste(title1,", p=",as.character(p1) ),col=col1,xlab=xlab1,ylab=ylab1, cex=cex1, cex.axis=cex.axis1,cex.lab=cex.lab1, cex.main=cex.main1)
}





# this function compute an asymptotic p-value for a given correlation (r) and sample size (n) 
if (exists("FisherTransformP") ) rm(FisherTransformP); 
FisherTransformP=function(r, n) {
#Z=sqrt(n-3) * 0.5*log((1+r)/(1-r)) 
#2*(1-pnorm(abs(Z) ))
# the following is implemented in the cor.test function
T=sqrt(n-2) * r/sqrt(1-r^2)
2*(1-pt(abs(T),n-2))
}


goodmargins1= c(0, 5, 2.3, 1)
