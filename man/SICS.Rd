\name{SICS}
\alias{SICS}
\title{Fit phylogeny-regularized GLM with phylogeny-constraint smoothness penalty}
\description{
    Fit a generalized linear model at grids of tuning parameter via penalized maximum likelihood.
    The regularization path is computed for a combination of sparse and smoothness penalty.
    The sparse penalty is MCP penalty and phylogeny-constraint smoothness penalty is derived from the phylogeny-induced correlation structure (C) based on patristic distance between OTU (D) in the form of C=exp(-2*pho*D). pho characterizes the evolutionary rate and performs phylogenetic grouping. Large pho groups OTUs into clusters at lower phylogenetic depth. The phylogeny-constraint smoothness penalty is in the form of beta^{T}*C^{-1}*beta. It has the advantage to encourage local smoothing by reducing the effects of distant OTUs. The estimator is called Sparse Inverse Correlation Srinkage Estimator (SICS).
}

\usage{
SICS<-function(X,Y,D,family=c("gaussian","binomial"),phos=c(2^(c(seq(-5,5,length=5)))),lambda2=c(0,2^(seq(-5, 5, length=5))),sparsity=0.9, nfolds=5)
}

\arguments{
  \item{X}{Input matrix; each row is an observation vector.}
  \item{Y}{Response vector. Quantitative for \code{family="gaussian"} or binary(0/1) for \code{family="binomial"}.}
   \item{D}{patristic distance between OTU.}
  \item{family}{Either "gaussian", "binomial", depending on the response.}
    \item{pho}{A user-specified sequence of \code{pho} values. pho characterizes the evolutionary rate and performs phylogenetic grouping. The default value is (0,2^(seq(-5, 5, length=5))).} 
   \item{lambda2}{A user-specified sequence of \code{lambda2} values.\code{lambda2} is the tuning parameter for smoothness penalty. The default value is (0,2^(seq(-5, 5, length=5))).}
   \item{nfolds}{folds to perform cross-validation. Default value is 5.}
}

\value{
	\item{beta}{A vector of fitted coefficients not include intercept}
}

\author{Li Chen <li.chen@auburn.edu>, Jun Chen <chen.jun2@mayo.edu>}
\examples{
	 
}
\keyword{models}
\keyword{regression}

