\name{MaternCorr}
\alias{Matern}
\alias{MaternCorr}
\alias{MaternCorr.default}
\alias{MaternCorr.dsCMatrix}
\alias{MaternCorr.dgCMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Matern correlation function and Matern formula term.
}
\description{
 The Matérn correlation function describes realizations of Gaussian spatial processes with different smoothnesses 
 (i.e. either smooth or rugged surfaces, controlled by the \eqn{\nu} parameter). It also includes a \eqn{\rho} scaling parameter and an optional 'nugget' parameter. A random effect specified in a model formula as \code{Matern(1|}\emph{<...>}\code{)} has pairwise correlations given by the Matérn function at the scaled Euclidean distance between coordinates specified in \emph{<...>}, using \dQuote{+} as separator (e.g., \code{Matern(1|longitude+latitude)}). The Matern family can be used in Euclidean spaces of any dimension; and also for correlations on a sphere (with maximum smoothness \code{nu=0.5}). 
 
A syntax of the form \code{Matern(1|longitude+latitude \%in\% grp)} can be used to specify a Matern random effect with independent realizations (but identical correlation parameters) for each level of the grouping variable \code{grp}. Alternatively, the \code{Matern(<T/F factor>|longitude+latitude)} may be used to specify Matern effects specific to individuals identified by the \code{<T/F factor>} (see Example with females and males). In that case distinct correlation parameters are fitted for each such Matern term.

When group-specific autocorrelated random effects are fitted, it may be wise to allow for different means for each group in the Intercept (a message will point this out if the fit results for Matern or Cauchy terms suggest so). 
 
By default, \code{fitme} and \code{corrHLfit} performs optimization over the \eqn{\rho} and \eqn{\nu} parameters. It is possible to estimate different scaling parameters for the different Euclidean dimensions: see examples in \code{\link{make_scaled_dist}}.
 
The \code{MaternCorr} function may be used to vizualise these correlations, using distances as input. 
}
\usage{
\method{MaternCorr}{default}(d, rho = 1, smoothness, nu = smoothness, Nugget = NULL)
# Matern(1|...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{A distance or a distance matrix.}
  \item{rho}{A scaling factor for distance. The 'range' considered in some 
            formulations is the reciprocal of this scaling factor}
  \item{smoothness}{The smoothness parameter, >0. \eqn{\nu=0.5} corresponds to the exponential correlation function, 
   and the limit function when \eqn{\mu} goes to \eqn{\infty} is the squared exponential function (as in a Gaussian).}
   
  \item{nu}{Same as smoothness}
  \item{Nugget}{(Following the jargon of Kriging) a parameter describing a discontinuous decrease in 
    correlation at zero distance. Correlation will always be 1 at \eqn{d=0}, and from which it immediately drops to 
    (1-Nugget)}
  \item{\dots}{Names of coordinates, using \dQuote{+} as separator (e.g., \code{Matern(1|longitude+latitude)}. The coordinates are numeric values found in the \code{data} data frame provided to the fitting function. No additional declaration of groups, factors, or other specific formatting is required.}
}
\details{
The correlation at distance \eqn{d>0} is 
  \deqn{(1-\textrm{Nugget}) \frac{(\rho d)^\nu  K_\nu(\rho d)}{2^{(\nu - 1)} \Gamma(\nu)}}{%
                         (1-Nugget)  2^(1-\nu) (\rho d)^\nu  K_\nu(\rho d) / \Gamma(\nu)} 
where
\eqn{K_\nu} is the \code{\link{besselK}} function of order \eqn{\nu}.

By default the Nugget is set to 0. See one of the examples on data set \code{\link{Loaloa}} 
for a fit including the estimation of the Nugget.
}
\value{Scalar/vector/matrix depending on input.}
\references{
Stein, M.L. (1999) Statistical Interpolation of Spatial Data: Some Theory for Kriging. Springer, New York.
}
\seealso{
See \code{\link{corMatern}} for an implementation of this correlation function as a \code{corSpatial} object for use with \code{lme} or \code{glmmPQL}.
}
\examples{
## See examples in help("HLCor"), help("Loaloa"),  help("make_scaled_dist"), etc.
## Matern correlations in 4-dimensional space:
set.seed(123)
randpts <- matrix(rnorm(20),nrow=5)
distMatrix <- as.matrix(proxy::dist(randpts))
MaternCorr(distMatrix,nu=2)

## Group-specific random effects: 
if (spaMM.getOption("example_maxtime")>1.6) {
  data(Leuca)
  subLeuca <- Leuca[c(1:10,79:88),] # subset of 10 females and 10 males, for faster examples.

  # Independent Matern random effect with different covariance parameters for each sex:
  fitme(fec_div ~ sex + Matern(female|x + y) + Matern(male|x + y),  data = subLeuca) 

  # Independent Matern random effect with the same covariance parameters for each sex:
  fitme(fec_div ~ sex + Matern(1|x+y \%in\% sex),data=subLeuca)

  # Matern models with random-effects distinct but correlated across sexes 
  # can also be fitted: see Matern examples in help("composite-ranef").   
} 
}
\keyword{ models }% use one of  RShowDoc("KEYWORDS")
\keyword{ spatial }% __ONLY ONE__ keyword per line
