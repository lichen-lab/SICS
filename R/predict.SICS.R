
predict.SICS<-function(beta,x,family){
  eta=as.numeric(cbind(1, x) %*% beta)
  yhat=eta
  if (family == 'binomial') {
    yhat=exp(eta) / (1 + exp(eta))
  }
  yhat
}




