inv.logit <- function (x) {
  exp(x) / (1 + exp(x))
}

#  lasso
lasso<-function(x,y,family='gaussian',nfolds=5){
  type.measure=ifelse(family=="binomial","auc","mse")
  cv=cv.glmnet(x,y,family=family,type.measure=type.measure, nfolds=nfolds)
  obj=glmnet(x,y,family=family,lambda=cv$lambda.min)
  beta=coef(obj)
  beta=as.vector(beta)
  beta
}


# MCP
MCP<-function(x,y,family='gaussian',nfolds=5){
  type.measure=ifelse(family=="binomial","auc","mse")
  cv=cv.ncvreg(x,y,family=family,penalty="MCP",type.measure=type.measure, nfolds=nfolds)	
  obj=ncvreg(x,y,family=family,penalty="MCP",lambda=cv$lambda.min)
  beta=coef(obj)
  beta=as.vector(beta)
  beta
}


# Elastic net
Enet<-function(x,y,family,nfolds=5,alphas=seq(0, 1, length=10)){
  type.measure=ifelse(family=="binomial","auc","mse")
  enets=lapply(alphas, function(alpha) { cv.glmnet(x, y, alpha=alpha, family=family,nfolds=nfolds) })
  enet=enets[[which.min(lapply(enets, function(enet) min(enet$cvm)))]]
  alpha=alphas[which.min(sapply(1:length(enets), function(x) min(enets[[x]]$cvm)))]
  obj=glmnet(x,y,family=family,alpha=alpha,lambda=enet$lambda.min)
  beta=coef(obj)
  beta=as.vector(beta)
  beta
}


# SLS
SLS<-function(x,y,family,nfolds,sparsity=0.9, D, pho=c(2^(c(seq(-5,5,length=5)))),lambda2=c(0,2^(seq(-5, 5, length=10)))){
  
  type.measure=ifelse(family=="binomial","auc","mse")
  betas=list()
  cvs=numeric()
  for (ph in pho){
    V=exp(-ph*D^2)
    VA=V[lower.tri(V)]
    SS<-quantile(VA, sparsity)
    V [V < SS]=0
    L= -V
    for (i in 1: ncol(D)){
      L[i,i]=sum(abs(V[i,]))-V[i,i]
    }
    cv=cv.glmgraph(x,y,L,family, penalty='MCP',lambda2=lambda2, type.measure=type.measure, nfolds=nfolds,trace=F)
    betas[[paste(ph)]]=cv$beta.min
    cvs[paste(ph)]=cv$cvmin
  }
  if(length(betas)==0){
    beta=NULL
  }else{
    cpho=intersect(names(betas),names(cvs))
    betas=betas[cpho]
    cvs=cvs[cpho]
    beta=betas[[which.min(cvs)]]
  }
  beta
}

eval0<-function(beta,z.te,y.te,family){
  
  p=ncol(z.te)
  pmse=R2=auc=NA
  F1=sens=spc=fdr=NA
  
  if(!is.null(beta)){
    eta=as.numeric(cbind(1, z.te) %*% beta)
    yhat=eta
    if (family == 'binomial') {
      yhat=inv.logit(eta)
    }
    pmse=mean((y.te-yhat)^2)
    R2=cor(y.te, yhat)
    if(family=="binomial"){
      #auc=calauc(y.te,yhat)
      auc=pROC::auc(y.te,yhat)
    }
  }
  res=t(data.frame(c(pmse=pmse,R2=R2,auc=auc)))
  rownames(res)=NULL
  res
}


pho=c(2^(seq(-5, 5, length=10)));lambda2=c(0,2^(seq(-5, 5, length=10)));nfolds=5;seed=1234
  

# repeated random sampling evaluation
eval<-function(x,y,D,family='gaussian',pho=c(2^(seq(-5, 5, length=10))),lambda2=c(0,2^(seq(-5, 5, length=10))),nrep=10,nfolds=5,seed=1234){
  
  set.seed(seed)
  
  method=c("SICS",'SLS(0)',"SLS(0.9)","Lasso","MCP",'Enet','RF')
  nmethod=length(method)
  crits=c('pmse','R2','auc')
  ncrit=length(crits)
  n=length(y)
  res=array(NA,c(nrep,nmethod,ncrit),dimnames=list(rep=1:nrep,method=method,crit=crits))

  for(irep in 1:nrep){
    message(irep)
    if(family=='gaussian'){
      cv.ind=ceiling(sample(1:n)/n*nfolds)
      x.tr=x[cv.ind!=1,]
      x.te=x[cv.ind==1,]
      y.tr=y[cv.ind!=1]
      y.te=y[cv.ind==1]
    }else if(family=='binomial'){
      cv.ind=ceiling(sample(1:n)/n*nfolds)
      ind1 <- which(y==1)
      ind0 <- which(y==0)
      n1 <- length(ind1)
      n0 <- length(ind0)
      cv.ind1 <- ceiling(sample(1:n1)/n1*nfolds)
      cv.ind0 <- ceiling(sample(1:n0)/n0*nfolds)
      cv.ind <- numeric(n)
      cv.ind[y==1] <- cv.ind1
      cv.ind[y==0] <- cv.ind0
      x.tr=x[cv.ind!=1,]
      x.te=x[cv.ind==1,]
      y.tr=y[cv.ind!=1]
      y.te=y[cv.ind==1]
    }
    
    beta.sics=SICS(X=x.tr,Y=y.tr,D=D,family=family,nfolds=nfolds,pho=pho,lambda2=lambda2)
    e.sics=eval0(beta.sics,x.te,y.te,family)
    res[irep,1,]=e.sics

    beta.sls=SLS(x.tr,y.tr,family,nfolds,sparsity=0, D, pho,lambda2)
    e.sls=eval0(beta.sls,x.te,y.te,family)
    res[irep,2,]=e.sls
    
    beta.sls=SLS(x.tr,y.tr,family,nfolds,sparsity=0.9, D, pho,lambda2)
    e.sls=eval0(beta.sls,x.te,y.te,family)
    res[irep,3,]=e.sls
    
    beta.ls=lasso(x.tr,y.tr,family,nfolds)
    e.ls=eval0(beta.ls,x.te,y.te,family)
    res[irep,4,]=e.ls
    
    beta.mcp=MCP(x.tr,y.tr,family,nfolds)
    e.mcp=eval0(beta.mcp,x.te,y.te,family)
    res[irep,5,]=e.mcp
    
    beta.enet=Enet(x.tr,y.tr,family,nfolds)
    e.enet=eval0(beta.enet,x.te,y.te,family)
    res[irep,6,]=e.enet
    
    
    e.rf=rep(NA,3)
    dd=data.frame(class=y.tr,x.tr)
    if(family=='binomial'){
      dd$class=as.factor(dd$class)
    }
    rf=randomForest(class~.,data=dd)
    colnames(x.te)=colnames(dd)[-1]
    if(family=='binomial'){
      yhat=predict(rf,x.te,type="prob")[,2]
      auc=pROC::auc(y.te,yhat)
      e.rf[3]=auc
    }else if(family=='gaussian'){
      yhat=predict(rf,x.te,type="response")
    }
    pmse=mean((y.te-yhat)^2)
    R2=cor(y.te, yhat)
    e.rf[1]=pmse
    e.rf[2]=R2
    res[irep,7,]=e.rf
    
  }
  
  res
}













