\name{HLCor}
\alias{HLCor}
\title{
                    Fits a (spatially) correlated mixed model, for given correlation parameters 
}
\description{
A convenient interface for \code{\link{HLfit}}, constructing the correlation matrix of random effects from the arguments, 
then estimating fixed effects and dispersion parameters using \code{HLfit}.
}
\usage{
HLCor(formula, data, family = gaussian(), ranPars = NULL, distMatrix,
      uniqueGeo = NULL, adjMatrix, corrMatrix, covStruct=NULL,
      method = "REML", verbose = c(trace=FALSE),
      control.dist = list(), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
   A \code{predictor}, i.e. a formula with attributes (see \code{\link{Predictor}}), or possibly simply a simple \code{formula} if an offset is not required.
  }
  \item{ranPars}{
  A list of values for correlation parameters (some of which are mandatory), and possibly also dispersion parameters 
  (optional, but passed to HLfit if present). See \code{\link{ranPars}} for further information.
}
  \item{data}{
   The data frame to be analyzed.
  }
  \item{family}{A \code{family} object describing the distribution of the response variable. See \code{\link{HLfit}} for further information.
  }
  \item{distMatrix}{
    A distance matrix between geographic locations, forwarded to \code{MaternCorr} or \code{CauchyCorr}. It overrides the (by default, Euclidean) distance matrix that would otherwise be deduced from the variables in a \code{Matern(.))} or \code{Cauchy(.)} term.   
}
  \item{uniqueGeo}{
    A matrix of non-redundant geographic locations. Such a matrix is typically constructed automatically from the data and the model formula, but otherwise could be useful if further the \code{rho} parameter is a vector with different values for different coordinates, in which case a scaled distance matrix has to be reconstructed from \code{uniqueGeo} and \code{rho}. 
}
  \item{adjMatrix}{
    An adjacency matrix, used if a random effect of the form\cr
    \code{y ~ adjacency(1|<location index>)} is present. See \code{\link{adjacency}} for further details.
}
  \item{corrMatrix}{
    A matrix \bold{C} used if a random effect term of the form \code{corrMatrix(1|<stuff>)} is present. This allows to analyze non-spatial model by giving for example a matrix of genetic correlations. Each row corresponds to levels of a variable <stuff>. The covariance matrix of the random effects for each level is then \eqn{\lambda}\bold{C}, where as usual \eqn{\lambda} denotes a variance factor for the random effects (if \bold{C} is a correlation matrix, then  \eqn{\lambda} is the variance, but other cases are possible). See \code{\link{corrMatrix}} for further details.
}
  \item{covStruct}{
    An interface for specifying correlation structures for different types of random effect (\code{corrMatrix} or \code{adjacency}). See \code{\link{covStruct}} for details.
  }
  \item{method}{Character: the fitting method to be used, such as \code{"ML"}, \code{"REML"} or \code{"PQL/L"}. \code{"REML"} is the default. Other possible values of \code{HLfit}'s \code{method} argument are handled.
}
  \item{verbose}{
    A vector of booleans. \code{trace} controls various diagnostic (possibly messy) messages about the iterations.
  }
  \item{control.dist}{
   A list of arguments that control the computation of the distance argument of the correlation functions. Possible elements are
   \describe{
   \item{rho.mapping}{ a set of indices controlling which elements of the \code{rho} scale vector scales which dimension(s) of the space in which (spatial) 
   correlation matrices of random effects are computed. See same argument in \code{\link{make_scaled_dist}} for details and examples.} 
   \item{dist.method}{\code{method} argument of \code{proxy::dist} function (by default, \code{"Euclidean"}, but see \code{\link{make_scaled_dist}} for other distances such as spherical ones.)}
   }
}
  \item{\dots}{
    Further parameters passed to \code{HLfit} or to \code{\link{mat_sqrt}}.
}
}
\details{
 For approximations of likelihood, see \code{\link{method}}.  For the possible structures of random effects, see \code{\link{random-effects}}, but note that \code{HLCor} cannot adjust parameters of correlation models (with the exception of conditional autoregressive ones). Any such parameter must be specified by the \code{ranPars} argument. More generally, the correlation matrix for random effects can be specified by various combinations of formula terms and other arguments  (see Examples):
    \describe{
    \item{Basic Matérn model}{\code{Matern(1|<...>)}, using the spatial coordinates in \code{<...>}. This will construct a correlation matrix according to the Matérn correlation function (see \code{\link{MaternCorr}});}
    \item{Basic Cauchy model}{\code{Cauchy(1|<...>)}, as for Matern (see \code{\link{CauchyCorr}});}
    \item{Same models with given distance matrix}{as provided by \code{distMatrix} (see Examples);}
    \item{Given correlation matrix}{\code{corrMatrix(1|<...>)} with \code{corrMatrix} argument. See \code{\link{corrMatrix}} for further details.}
    \item{CAR model with given adjacency matrix}{\code{adjacency(1|<...>)} with \code{adjMatrix}. See \code{\link{adjacency}} for further details;}
    \item{AR1 model}{\code{AR1(1|<...>)} See \code{\link{AR1}} for further details.} 
  }
}

\value{
  The return value of an \code{HLfit} call, with the following additional attributes: 
  \item{HLCorcall}{the HLCor call}  
  \item{info.uniqueGeo}{Unique geographic locations.}  
}

\seealso{
\code{\link{autoregressive}} for additional examples, \code{\link{MaternCorr}}, \code{\link{HLfit}}, and \code{\link{corrHLfit}}
}
\examples{
# Example with an adjacency matrix (autoregressive model):
# see 'adjacency' documentation page

#### Matern correlation using only the Matern() syntax
data("blackcap")
(fitM <- HLCor(migStatus ~ means+ Matern(1|latitude+longitude),data=blackcap,
      method="ML", ranPars=list(nu=0.6285603,rho=0.0544659)))

#### Using the 'distMatrix' argument
data("blackcap")
#
# Build distance matrix (here equivalent to the default one for a Matern() term)
MLdistMat <- as.matrix(proxy::dist(blackcap[,c("latitude","longitude")]))
#
(fitD <- HLCor(migStatus ~ means+ Matern(1|latitude+longitude),data=blackcap,
      distMatrix=MLdistMat, method="ML", ranPars=list(nu=0.6285603,rho=0.0544659)))
# : result here must be equivalent to the one without the distMatrix.
diff(c(logLik(fitM),logLik(fitD)))
}
\keyword{ model }
