\name{SICS-package}
\alias{SICS-package}
\docType{package}
\title{Fit a GLM with a combination of sparse and phylogeny-induced smoothness penalty}
\description{
    Fit a generalized linear model at grids of tuning parameter via penalized maximum likelihood.
    The regularization path is computed for a combination of sparse and smoothness penalty.
    The sparse penalty is MCP penalty and the smoothness penalty has two options. The phylogeny-induced correlation structure (C) is derived from patristic distance between OTU (D) in the form of C=exp(-2*pho*D). pho characterizes the evolutionary rate and performs phylogenetic grouping. Large pho groups OTUs into clusters at lower phylogenetic depth. If the smoothness penalty is the phylogeny-constraint smoothness penalty, it is in the form of beta^{T}*C^{-1}*beta. It has the advantage to encourage local smoothing by reducing the effects of distant OTUs. The estimator is called Sparse Laplacian Srinkage Estimator (SLS). If the smoothness penalty is the
   Laplacian penalty, it is in the form of beta^{T}*L*beta. By sparsifying C, we obtain Laplacian matrix (L) with different level of sparsity. The estimator is called Sparse Inverse Correlation Srinkage Estimator (SICS).
}
\details{
  \tabular{ll}{
    Package: \tab SICS\cr
    Type: \tab Package\cr
    Version: \tab 1.0-0\cr
    Date: \tab 2018-08-27\cr
    License: \tab GPL-2\cr}
  The algorithm accepts a design matrix \code{X}, a vector of responses \code{Y} and a distance matrix \code{D}.
  It consists of the following main functions \cr
  \code{glmgraph}\cr
  \code{cv.glmgraph}\cr
   \code{SLS}\cr
  \code{SICS}\cr
  }

\author{Li Chen <li.chen@auburn.edu>, Jun Chen <jun.chen2@mayo.edu>}
\keyword{models}
\keyword{regression}
\keyword{package}
\examples{
 set.seed(1234)
 library(SICS)
  
}









