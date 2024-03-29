\name{make_scaled_dist}
\alias{make_scaled_dist}
\alias{rho.mapping}
\alias{Earth}
\alias{EarthChord}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Scaled distances between unique locations}
\description{
  This function computes scaled distances from whichever relevant argument it can use (see Details). The result can directly by used as input for computation of the Matérn correlation matrix. It is usually called internally by HLCor, so that users may ignore it, except if they wish to control the distance used through \code{control.dist$method}, or the parametrization of the scaling through \code{control.dist$rho.mapping}. \code{control.dist$method} provide access to the distances implemented in the \code{proxy} package, as well as to \code{"EarthChord"} and \code{"Earth"} methods defined in \code{spaMM} (see Details).
}
\usage{
make_scaled_dist(uniqueGeo, uniqueGeo2=NULL, distMatrix, rho, 
                 rho.mapping=seq_len(length(rho)), 
                 dist.method="Euclidean",
                 return_matrix=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{uniqueGeo}{
     A matrix of geographical coordinates (e.g. 2 columns for latitude and longitude), without replicates of the same location.  
}
  \item{uniqueGeo2}{
     NULL, or a second matrix of geographical coordinates, without replicates of the same location.  
     If NULL, scaled distances among \code{uniqueGeo} locations are computed. Otherwise, scaled distances between  locations in the two input matrices are computed. 
}
  \item{distMatrix}{
    A distance matrix.
}
  \item{rho}{
   A scalar or vector of positive values. Scaled distance is computed as \code{<distances in each coordinate> * rho}, unless a non-trivial
   \code{rho.mapping} is used.
}
  \item{rho.mapping}{
     A set of indices controlling which elements of the \code{rho} scale vector scales which dimension(s) of the space in which (spatial) 
   correlation matrices of random effects are computed. Scaled distance is generally computed as \code{<distances in each coordinate> * rho[rho.mapping]}. As shown in the Example, if one wishes to combine isotropic geographical distance and some environmental distance, the coordinates being latitude, longitude and one environmental variable, the scaled distance may be computed
  as (say) \code{(lat,long,env) *rho[c(1,1,2)]} so that the same scaling \code{rho[1]} applies for both geographical coordinates. In this case, \code{rho} should have length 2 and \code{rho.mapping} should be \code{c(1,1,2)}.
  }
  \item{dist.method}{\code{method} argument of \code{proxy::dist} function (by default, \code{"Euclidean"}, but other distances are possible (see Details).} 
  \item{return_matrix}{Whether to return a \code{matrix} rather than a \code{proxy::dist} or \code{proxy::crossdist} object.}
}
\details{
 The function uses the \code{distMatrix} argument if provided, in which case rho must be a scalar. Vectorial \code{rho} (i.e., different scaling of different dimensions) is feasible only by providing \code{uniqueGeo}.
 
 The \code{dist.method} argument gives access to distances implemented in the \code{proxy} package, or to user-defined ones that are made accessible to \code{proxy} through its database. Of special interest for spatial analyses are distances computed from longitude and latitude (\code{proxy} implements \code{"Geodesic"} and \code{"Chord"} distances but they do not use such coordinates: instead, they use Euclidean distance for 2D computations, i.e. Euclidean distance between points on a circle rather than on a sphere). spaMM implements two such distances: \code{"Earth"} and \code{"EarthChord"}, using longitude and latitude inputs \bold{in that order} (see Examples). The \code{"EarthChord"} distance is the 3D Euclidean distance \dQuote{through Earth}. The \code{"Earth"} distance is also known as the orthodromic or great-circle distance, on the Earth surface. Both distances return values in km and are based on approximating the Earth by a sphere of radius 6371.009 km.  
}
\value{
A matrix or \code{\link[proxy]{dist}} object. If there are two input matrices, rows of the return value correspond to rows of the first matrix.
}
\examples{
data("blackcap")
## a biologically not very meaningful, but syntactically correct example of rho.mapping
fitme(migStatus ~ 1 + Matern(1|longitude+latitude+means),
      data=blackcap, fixed=list(nu=0.5,phi=1e-6),
      init=list(rho=c(1,1)), control.dist=list(rho.mapping=c(1,1,2)))

## Using orthodromic distances: 
# order of variables in Matern(.|longitude+latitude) matters; 
# Matern(1|latitude+longitude) should cause a warning
fitme(migStatus ~ 1 + Matern(1|longitude+latitude),data=blackcap,
      method="ML", fixed=list(nu=0.5,phi=1e-6),
      control.dist=list(dist.method="Earth"))
}
\keyword{ model }
