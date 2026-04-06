\name{optimBounds}
\alias{optimBounds}
\alias{lower}
\alias{upper}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Optimization bounds
}
\description{
Higher-level fitting functions (\code{fitme}, its derived function \code{fitmv} and \code{pois4mlogit}, and the older function \code{corrHLfit}), perform box-constrained numerical optimization of the likelihood objective function. The constraints can be specified in calls of these functions by the \code{lower} and \code{upper} arguments. 

Specifying upper bounds for variances of random effects, and lower and upper bounds for the correlation parameters of random-coefficient covariance matrices, may be useful to prevent divergence of their estimates. Such divergence may otherwise occur notably in REML or ML fits of binary-response models, and more broadly of count-response data with low count values.

The experimental \code{optimBounds} function may be used to extract information from a fit object about the lower and upper bounds used in optimization.
}
\usage{
optimBounds(x, transf, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Fit object produced by \code{fitme} or related higher-level fitting functions.
}
  \item{transf}{
  boolean; whether to report transformed values of (co-)variance parameters or not. Random-coefficients parameter bounds are always returned in transformed space, for reasons given in Details. 
}
  \item{\ldots}{
  For development purposes.
}
}
\details{
Optimization of random-coefficient parameters is performed in a transformed parameter space. Box constraints on the elements of the \code{ranCoefs} vectors (i.e., variances and correlations), do not translate into box constraints in the transformed space, and vice versa. For this reason, simple transformation of lower or upper bounds from \code{ranCoefs} space to transformed space is generally inadequate, and \code{optimBounds} will not provide random-coefficient bounds untransformed as \code{ranCoefs}, even when \code{transf=FALSE}. It will instead keep any transformed values \code{trRanCoefs} in the result, with a message. Such \code{trRanCoefs} do not take into account optional box constraints declared by the user as \code{ranCoefs} elements of the \code{lower} or \code{upper} arguments of the fitting function that returned the fit object. Yet, since version 4.6.36, the fit object does satisfy such constraints correctly.   
}
\value{
A structured list, with elements \code{lower} and \code{upper} and themselves being lists.   
}
\examples{
# Example from VarCorr() documentation in 'nlme' package:
data("Orthodont",package = "nlme")
sp1 <- fitme(distance ~ age+(age|Subject), data = Orthodont, method="REML")
optimBounds(sp1, transf=TRUE) # showing transformed values 'trRanCoefs'

# Constrained fit (here with clear effect on first variance):
sp2 <- fitme(distance ~ age+(age|Subject), data = Orthodont, method="REML",
             upper=list(ranCoefs=list("1"=c(1,0.999,1))))
optimBounds(sp2, transf=TRUE) # showing transformed values 'trRanCoefs' 

# Specifying bounds on transformed parameters is possible:
sp3 <- fitme(distance ~ age+(age|Subject), data = Orthodont, method="REML",
             upper=list(trRanCoefs=list("1"=c(1,0.999,0.0447))))
optimBounds(sp3, transf=TRUE) 

}
