SICS<-function(X,Y,D,family=c("gaussian","binomial"),pho=c(2^(c(seq(-5,5,length=5)))),lambda2=c(0,2^(seq(-5, 5, length=10))), nfolds=5){
  
  if(missing(X) | missing(Y) | missing (D)) stop("x, y, D are required!")
  type.measure=ifelse(family=="binomial","auc","mse")
  betas=list()
  cvs=numeric()
  for (p in pho){
    V=exp((-1)*p *D^2)
    L=as.matrix(nearPD(ginv(V))$mat)
    VL=matrix(0, nrow=ncol(D), ncol=ncol(D))
    for (i in 1: ncol(D)){
      VL[i,i]=(sum(abs(L[i,])))^(-1/2)
    }
    L=VL %*% L %*% VL
    cv=cv.glmgraph(X,Y,L,family, penalty='MCP',lambda2=lambda2, type.measure=type.measure, nfolds=nfolds,trace=F)
    betas[[paste(p)]]=cv$beta.min
    cvs[paste(p)]=cv$cvmin
  }
  
  if(length(betas)==0){
    beta=NULL
  }else{
    cp=intersect(names(betas),names(cvs))
    betas=betas[cp]
    cvs=cvs[cp]
    beta=betas[[which.min(cvs)]]
  }

  structure(beta, class="SICS")
  
}








