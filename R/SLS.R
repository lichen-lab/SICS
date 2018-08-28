SLS<-function(X,Y,D,family=c("gaussian","binomial"),pho=c(2^(c(seq(-5,5,length=5)))),lambda2=c(0,2^(seq(-5, 5, length=10))), nfolds=5,sparsity=0.9){
  
  if(missing(X) | missing(Y) | missing (D)) stop("x, y, D are required")
  type.measure=ifelse(family=="binomial","auc","mse")
  betas=list()
  cvs=numeric()
  for (p in pho){
    V=exp((-1)*p *D^2)
    VA=V[lower.tri(V)]
    SS<-quantile(VA, sparsity)
    V [V < SS]=0
    L= -V
    for (i in 1: ncol(D)){
      L[i,i]=sum(abs(V[i,]))-V[i,i]
    }
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
  structure(beta, class="SLS")
  
}














