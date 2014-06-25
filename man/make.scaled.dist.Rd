\name{make.scaled.dist}
\alias{make.scaled.dist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Scaled Euclidian distances between unique locations}
\description{
  This fonction computes scaled Euclidian distances from whichever relevant argument it can use (see Details). The result can directly by used as input for computation of the Matérn correlation matrix. It is usually called internally by HLCor, so that users may ignore it, except if they wish to control 
  the parametrizaion of the scaling through the \code{rho.mapping} argument.
}
\usage{
make.scaled.dist(distMatrix,uniqueGeo,rho,rho.mapping=seq_len(length(rho)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{distMatrix}{
    A distance matrix.
}
  \item{uniqueGeo}{
     A matrix of geographical coordinates (e.g. 2 columns for latitude and longitude), without replicates o the same location.  
}
  \item{rho}{
   A scalar or vector of positive values. Scaled distance is computed as \code{<distances in each coordinate> * rho}, unless a non-trivial
   \code{rho.mapping} is used.
}
  \item{rho.mapping}{
     A set of indices controlling which elements of the \code{rho} scale vector scales which dimension(s) of the space in which (spatial) 
   correlation matrices of random effects are computed. Scaled distance is generally computed as \code{<distances in each coordinate> * rho[rho.mapping]}. As shown inthe Example, if one wishes to combine isotropic geographical distance and some environmental distance, the coordinates being latitude, longitude and one environmental variable, the scaled distance may be computed
  as (say) \code{(lat,long,env) *rho[c(1,1,2)]} so that the same scaling \code{rho[1]} applies for both geographical coordinates. In this case, \code{rho} should have length 2 and \code{rho.mapping} should be \code{c(1,1,2)}.
  }
}
\details{
 The function uses the \code{distMatrix} argument if provided, in which case rho must be a scalar. Vectorial \code{rho} (i.e., different scaling of different dimensions) is feasible only by providing \code{uniqueGeo}.
}
\value{
A matrix or \code{\link[proxy]{dist}} object.
}
\examples{
## a biologically not very meaningful, but syntactically correct example of rho.mapping
data(blackcap)
corrHLfit(migStatus ~ 1 + Matern(1|latitude+longitude+means),data=blackcap,
          objective="p_v",HLmethod="ML",ranFix=list(nu=0.5),init.corrHLfit=list(rho=c(1,1)),
          rho.mapping=c(1,1,2))
}
\keyword{ model }
