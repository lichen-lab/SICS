library(Rcpp)
library(dirmult)
library(ape)
library(ade4)
library(cluster)
library(MASS)

library(SICS)
library(randomForest)
library(glmgraph)
library(glmnet)
library(ncvreg)
library(vegan)
library(GMPR)



load('data_caff.RData')
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
save(data_caff_norm,file='data_caff_norm.RData')


### perform CV evaluation






















codepath='/scratch/bioinfo2/lichen/AOAS/realdata/codes'
datapath='/scratch/bioinfo2/lichen/AOAS/realdata/data'

source(file.path(codepath,'func.R'))






