\name{SLS}
\alias{SLS}
\title{Fit phylogeny-regularized GLM with Laplacian smoothness penalty}
\description{
    Fit a generalized linear model at grids of tuning parameter via penalized maximum likelihood.
    The regularization path is computed for a combination of sparse and smoothness penalty.
    The sparse penalty is MCP penalty and the smoothness penalty is Laplacian penalty derived from the phylogeny-induced correlation structure (C) based on patristic distance between OTU (D) in the form of C=exp(-2*pho*D). pho characterizes the evolutionary rate and performs phylogenetic grouping. Large pho groups OTUs into clusters at lower phylogenetic depth. By sparsifying C, we obtain Laplacian matrix (L) with different level of sparsity. The Laplacian penalty is in the form of beta^{T}*L*beta.
    Fits linear, logistic regression models. The estimator is called Sparse Laplacian Srinkage Estimator (SLS).
}

\usage{
SLS(X,Y,D,family=c("gaussian","binomial"),phos=c(2^(c(seq(-5,5,length=5)))),lambda2=c(0,2^(seq(-5, 5, length=5))),sparsity=0.9, nfolds=5)
}

\arguments{
  \item{X}{Input matrix; each row is an observation vector.}
  \item{Y}{Response vector. Quantitative for \code{family="gaussian"} or binary(0/1) for \code{family="binomial"}.}
   \item{D}{patristic distance between OTU.}
  \item{family}{Either "gaussian", "binomial", depending on the response.}
    \item{pho}{A user-specified sequence of \code{pho} values. pho characterizes the evolutionary rate and performs phylogenetic grouping. The default value is (0,2^(seq(-5, 5, length=5))).} 
   \item{lambda2}{A user-specified sequence of \code{lambda2} values.\code{lambda2} is the tuning parameter for smoothness penalty. The default value is (0,2^(seq(-5, 5, length=5))).}
   \item{sparsity}{The sparsity level to convert \code{C} to \code{L}. It is between 0 and 1. }
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

