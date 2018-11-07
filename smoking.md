
**1. Introduction of dataset**
-------------------------------
The smoking data was from a study of the smoking effect on the human upper respiratory tract micro-
384 biome. We aimed to predict the smoking status based on the microbiome profile. All
385 the data processing steps were carried out as described in the previous example. After preprocessing, the
386 final dataset consisted of 32 non-smokers and 28 smokers with 174 OTUs.

**2. Data processing**
-------------------
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


# load caffeine data
data('data_smoking')
otu.tab=data_smoking$otu.tab
smoking=data_smoking$smoking
tree=data_smoking$tree


# step1: filter OTU
n=nrow(otu.tab)
threshold=0.9
zero.per = colSums(otu.tab==0)/n
nzero.idx = which(zero.per<threshold)
zero.idx=which(zero.per>=threshold)
otu.ids1=colnames(otu.tab)[nzero.idx]
otu.median=apply(otu.tab,2,function(x) median(x[x!=0]))
otu.median[is.na(otu.median)]=0
otu.ids2=colnames(otu.tab)[otu.median>0]
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
y=as.numeric(smoking)-1
data_smoking_norm=NULL
data_smoking_norm$x=x
data_smoking_norm$y=y
data_smoking_norm$D=D

# step5: perform random-sampling evaluation
x=data_caff_norm$x
y=data_caff_norm$y
D=data_caff_norm$D
```

**3. Method comparison**
-------------------
To have an objective evaluation of the prediction performance, the dataset was randomly divided fifty
times into five folds each time, among which four folds were used for training and the remaining one for
testing. In the training set, tuning parameter selection was based on CV as in the simulation. R2 and PMSE
were used as metrics for prediction performance based on the testing set. SLS, Lasso, MCP, Elastic Net and
Random Forest are compared prediction methods.
```
# Download library.R from github website, which contains functions of competing prediction methods and prediction assessment functions
source('library.R')
nrep=1 # dataset was randomly divided nrep, default is 50
res=eval(x,y,D,family='binomial',nrep=nrep,nfolds=5,seed=1234)
res
```
