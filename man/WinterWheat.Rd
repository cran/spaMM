\name{WinterWheat}
\alias{WinterWheat}
\alias{GxE}
\title{
Example of yield stability analysis
}
\description{
Translation of an example that may be found at\cr
\url{https://www.r-bloggers.com/2019/06/genotype-experiments-fitting-a-stability-variance-model-with-r/},
based on yield of eight durum wheat genotypes over seven years, following a randomised block design with three replicates.
A genotype-in-year random effect is used to quantify genotype-by-environment interactions. 
In the first fit (\code{constvar}), the variance of this random effect is constant over genotypes. 
In the second fit (\code{varvar}), different variances are fitted for the distinct genotypes, to assess the relative stability of yield of the different genotypes over environments. This second model can be fitted as a 
constrained random-coefficient model, where the constraint describes a diagonal covariance matrix for the random coefficients.

This example uses the fact that the argument \code{fixed=list(ranCoefs=<...>)} can be used to fit a 
covariance matrix with an arbitrary set of constrained elements. Only elements 
left as 'NA' (here the diagonal elements of the matrix) are fitted.  

}
\examples{

if (spaMM.getOption("example_maxtime")>1.5 &&
      requireNamespace("agridat", quietly = TRUE)) {

data("onofri.winterwheat", package="agridat")

(constvar <- fitme(
  yield ~ gen + (1|year) + (1|block \%in\% year)+(1|gen \%in\% year), 
  data=onofri.winterwheat, method="REML"))
  
# Diagonal matrix of NA's, represented as vector for its lower triangle:
ranCoefs_for_diag <- function(nlevels) { 
 ## Conceptual version
 # diagmat <- matrix(NA, ncol=nlevels,nrow=nlevels)
 # diagmat[lower.tri(diagmat,diag=FALSE)] <- 0
 # return(diagmat[lower.tri(diagmat,diag=TRUE)])
 ## which amount to:
 vec <- rep(0,nlevels*(nlevels+1L)/2L)
 vec[cumsum(c(1L,rev(seq(nlevels-1L)+1L)))] <- NA
 vec
} 

(varvar <- fitme(
  yield ~ gen + (1|year) + (1|block \%in\% year)+(0+gen|gen \%in\% year), method="REML", 
  data=onofri.winterwheat, fixed=list(ranCoefs=list("3"=ranCoefs_for_diag(8L)))))

}

}
