\name{MaternCorr}
\alias{Matern}
\alias{MaternCorr}
\alias{MaternCorr.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Matern correlation function and Matern formula term.
}
\description{
 The Matérn correlation function describes realizations of Gaussian spatial processes with different smoothnesses 
 (i.e. either smooth or rugged surfaces). It also includes a scaling and a 'nugget' parameter. It can be invoked in two ways. First, the \code{MaternCorr} function evaluates these correlations, using distances as input. Second, a term of the form \code{Matern(1|}\emph{<...>}\code{)} in a \code{formula} specifies a random effect with Matérn correlation function, using coordinates found in a data frame as input. In the latter case, the correlations between realizations of the random effect for any two observations in the data will be the value of the Matérn function at the scaled Euclidean distance between coordinates specified in \emph{<...>}, using \dQuote{+} as separator (e.g., \code{Matern(1|latitude + longitude)}). 
}
\usage{
\method{MaternCorr}{default}(d, rho = 1, smoothness, nu = smoothness, Nugget = 0L)
# Matern(1|...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{A distance, typically an Euclidean distance}
  \item{rho}{A scaling factor for distance. The 'range' considered in some 
            formulations is the reciprocal of this scaling factor}
  \item{smoothness}{The smoothness parameter, >0. \eqn{\nu=0.5} corresponds to the exponential correlation function, 
   and the limit function when \eqn{\mu} goes to \eqn{\infty} is the squared exponential function (as in a Gaussian).}
   
  \item{nu}{Same as smoothness}
  \item{Nugget}{(Following the jargon of Kriging) a parameter describing a discontinuous decrease in 
    correlation at zero distance. Correlation will always be 1 at \eqn{d=0}, and from which it immediately drops to 
    (1-Nugget)}
  \item{...}{Names of coordinates, using \dQuote{+} as separator (e.g., \code{Matern(1|latitude + longitude)}}
}
\details{
The correlation at distance \eqn{d>0} is 
  \deqn{(1-\textrm{Nugget}) \frac{(\rho d)^\nu  K_\nu(\rho d)}{2^{(\nu - 1)} \Gamma(\nu)}}{%
                         (1-Nugget)  2^(1-\nu) (\rho d)^\nu  K_\nu(\rho d) / \Gamma(\nu)} 
where
\eqn{K_\nu} is the \code{\link{besselK}} function of order \eqn{\nu}.
}
\value{Scalar/vector/matrix depending on input.}
\references{
Stein, M.L. (1999) Statistical Interpolation of Spatial Data: Some Theory for Kriging. Springer, New York.
}
\seealso{
See \code{\link{corMatern}} for an implementation of this correlation function as a \code{corSpatial} object for use with \code{lme} or \code{glmmPQL}.

By default the Nugget is set to 0. See one of the examples on data set \code{\link{Loaloa}} 
for a fit including the estimation of the Nugget.
}
\examples{
## See examples in help("spaMM"), help("HLCor"), help("Loaloa"), etc.
## The Matérn function can be used in Euclidean spaces of any dimension:
set.seed(123)
randpts <- matrix(rnorm(20),nrow=5)
distMatrix <- as.matrix(proxy::dist(randpts))
MaternCorr(distMatrix,nu=2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
