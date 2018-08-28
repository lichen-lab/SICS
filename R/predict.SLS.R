
predict.SLS<-function(beta,X,family){
  eta=as.numeric(cbind(1, X) %*% beta)
  yhat=eta
  if (family == 'binomial') {
    yhat=exp(eta) / (1 + exp(eta))
  }
  yhat
}

