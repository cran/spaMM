\name{Matern.corr}
\alias{Matern.corr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Matern correlation function.
}
\description{
 A family of correlation function describing realizations of Gaussian spatial processes with different smoothnesses 
 (i.e. either smooth or rugged surfaces). Also includes a scaling and a 'nugget' parameter. 
}
\usage{
Matern.corr(d, rho = 1, smoothness, nu = smoothness, Nugget = 0L)
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
By default the Nugget is set to 0. See one of the examples on data set \code{\link{Loaloa}} 
for a fit including the estimation of the Nugget.
}
\examples{
## The Matérn function can be used in Euclidean spaces of any dimension:
set.seed(123)
randpts <- matrix(rnorm(20),nrow=5)
distMatrix <- as.matrix(proxy::dist(randpts))
Matern.corr(distMatrix,nu=2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
