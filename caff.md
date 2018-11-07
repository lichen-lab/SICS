The dataset was downloaded from Qiita (https://qiita.ucsd.edu/) with study ID 1011, which consists of 98 samples and 6674 OTUs. We selected the caffeine intake as the outcome of interest since caffeine intake was found to have a significant impact on the gut microbiota. We aimed to predict the caffeine intake based on the OTU abundances. Before applying the prediction methods, we implemented a series of preprocessing steps designed in to make the microbiome data more amenable to predictive modeling. First, we removed outlier samples based on an outlier index defined on Bray-Curtis distance and removed rare OTUs with prevalence less than 10% to reduce the dimensionality of OTUs, leaving 98 samples and 499 OTUs. Second, we normalized OTU raw read counts using GMPR followed by a replacement of outlier counts using winsorization at 97% quantile. Third, we transformed the normalized OTU abundance data using square-root transformation to reduce the influence of highly abundant observation. Finally, we applied quantile transformation to the caffeine intake to make it approximately normally distributed.



```
library(ape)
library(ade4)
library(cluster)
library(randomForest)
library(glmgraph)
library(glmnet)
library(ncvreg)
library(SICS)
library(GMPR)


# contains functions of competing prediction methods and prediction assessment functions
source("library.R") 
# load caffeine data
data('data_caff')
otu.tab=data_caff$otu.tab
caff=data_caff$caff
tree=data_caff$tree


# step1: filter OTU

n=nrow(otu.tab)
threshold=0.9
zero.per = colSums(otu.tab==0)/n
nzero.idx = which(zero.per<threshold)
zero.idx=which(zero.per>=threshold)
otu.ids1=colnames(otu.tab)[nzero.idx]
otu.median=apply(otu.tab,2,function(x) median(x[x!=0]))
otu.median[is.na(otu.median)]=0
otu.ids2=colnames(otu.tab)[otu.median>1]
otu.ids=intersect(otu.ids1,otu.ids2)
tree.tips=tree$tip.label
common.tips=intersect(tree.tips,otu.ids); length(common.tips)
tree=drop.tip(tree, setdiff(tree.tips, common.tips))
D=cophenetic(tree)
otu.ids=common.tips
otu.tab=otu.tab[,common.tips]

# step2: normalization
gmpr.size.factor=GMPR(otu.tab)
summary(gmpr.size.factor)
otu.tab.norm= otu.tab / gmpr.size.factor

# step3: remove outlier
q97=apply(otu.tab.norm,2,function(x) quantile(x,0.97))
otu.tab.win=apply(otu.tab.norm,2,function(x) {x[x>quantile(x,0.97)]=quantile(x,0.97); x} )
x=otu.tab.win

# step4: transformation
x[x!=0]=sqrt(x[x!=0])
y=qqnorm(caff)$x

data_caff_norm=NULL
data_caff_norm$x=x
data_caff_norm$y=y
data_caff_norm$D=D

# step5: perform random-sampling evaluation

x=data_caff_norm$x
y=data_caff_norm$y
D=data_caff_norm$D

res=eval(x,y,D,family='gaussian',nrep=2,nfolds=5,seed=1234)
res

```
