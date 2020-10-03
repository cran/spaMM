\name{CauchyCorr}
\alias{CauchyCorr}
\alias{Cauchy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cauchy correlation function and Cauchy formula term
}
\description{
 The Cauchy family of correlation functions is useful to describe spatial processes with power-law decrease of correlation at long distance. It is valid for Euclidean distances in spaces of any dimension, and for great-circle distances on spheres of any dimension. It has a scale parameter (\code{rho}, as in the Matérn correlation function), a \code{shape} (or \dQuote{smoothness}, Gneiting 2013) parameter, and a \code{long}-memory \code{dep}endence (or, more abstractly, \dQuote{shape}; Gneiting 2013) parameter (Gneiting and Schlater 2004). The present implementation also accepts a \code{Nugget} parameter. The family can be invoked in two ways. First, the \code{CauchyCorr} function evaluates correlations, using distances as input. Second, a term of the form \code{Cauchy(1|}\emph{<...>}\code{)} in a \code{formula} specifies a random effect with Cauchy correlation function, using coordinates found in a data frame as input. In the latter case, the correlations between realizations of the random effect for any two observations in the data will be the value of the Cauchy function at the scaled distance between coordinates specified in \emph{<...>}, using \dQuote{+} as separator (e.g., \code{Cauchy(1|latitude + longitude)}). A syntax of the form \code{Cauchy(1|latitude + longitude \%in\% grp)} can be used to specify a Cauchy random effect with independent realizations for each level of the grouping variable \code{grp}. 
}
\usage{
CauchyCorr(d, rho=1, shape, longdep, Nugget=NULL)
# Cauchy(1|...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{Euclidean or great-circle distance}
  \item{rho}{The scaling factor for distance, a real >0.}
  \item{shape}{The shape (smoothness) parameter, a real 0<.<=2 for Euclidean distances and 0<.<=1 for great-circle distances. Smoothness increases, and fractal dimension decreases, with increasing shape (the fractal dimension of realizations in spaces of dimension \eqn{d} being \eqn{d}+1-shape/2).}
  \item{longdep}{The long-memory dependence parameter, a real >0. It gives the exponent of the asymptotic decrease of correlation with distance: the \bold{smaller} \code{longdep} is, the longer the dependence.}
  \item{Nugget}{(Following the jargon of Kriging) a parameter describing a discontinuous decrease in 
    correlation at zero distance. Correlation will always be 1 at \eqn{d=0}, and from which it immediately drops to 
    (1-Nugget). Defaults to zero.}
  \item{\dots}{Names of coordinates, using \dQuote{+} as separator (e.g., \code{Matern(1|latitude + longitude)}}
}
\details{
The correlation at distance \eqn{d>0} is 
  \deqn{(1-\textrm{Nugget}) (1+(\rho d)^textrm{shape})^(-textrm{longdep/shape})}{%
                         (1-Nugget)  (1+(\rho d)^shape)^(-longdep/shape)} 
}
\value{Scalar/vector/matrix depending on input.}
\references{
Gneiting, T. and Schlater M. (2004) Stochastic models that separate fractal dimension and the Hurst effect. SIAM Rev.
46: 269–282.

Gneiting T. (2013) Strictly and non-strictly positive definite functions on spheres. Bernoulli 19: 1327-1349.
}
\examples{
data("blackcap")
fitme(migStatus ~ means+ Cauchy(1|latitude+longitude),data=blackcap,
      fixed=list(longdep=0.5,shape=0.5,rho=0.05))
## The Cauchy family can be used in Euclidean spaces of any dimension:
set.seed(123)
randpts <- matrix(rnorm(20),nrow=5)
distMatrix <- as.matrix(proxy::dist(randpts))
CauchyCorr(distMatrix,rho=0.1,shape=1,longdep=10)

# See ?MaternCorr for examples of syntaxes for group-specific random effects,
#  also handled by Cauchy().
}
\keyword{ models }% use one of  RShowDoc("KEYWORDS")
\keyword{ spatial }% __ONLY ONE__ keyword per line
