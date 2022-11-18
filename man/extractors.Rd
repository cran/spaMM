\name{extractors}
\alias{model.frame.HLfit} % undocumented as long as it is not stable (this has limitations)
\alias{model.matrix.HLfit} % not used in spaMM code (?), but possibly by other post-fit packages; and undocumented
\alias{extractors}
\alias{getDistMat}
\alias{logLik}
\alias{logLik.HLfit}
\alias{fitted}
\alias{fitted.HLfit}
\alias{fixef}
\alias{fixef.HLfit}
\alias{formula}
\alias{formula.HLfit}
\alias{family}
\alias{family.HLfit}
\alias{terms}
\alias{terms.HLfit}
\alias{nobs}
\alias{nobs.HLfit}
\alias{ranef}
\alias{ranef.HLfit}
\alias{print.ranef}
\alias{deviance}
\alias{deviance.HLfit}
\alias{response}
\alias{dev_resids}
\alias{df.residual}
\alias{df.residual.HLfit}
\alias{weights}
\alias{weights.HLfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Functions to extract various components of a fit.
}
\description{
 Most extractors are methods of generic functions defined in base R (see Usage), for which the base documentation may be useful.
 
\code{formula} extracts the model formula.\cr
\code{family} extracts the response family.\cr
\code{terms} extracts the formula, with attributes describing the \bold{fixed-effect} terms.\cr 
\code{nobs} returns the length of the response vector.\cr
\code{logLik} extracts the log-likelihood (exact or approximated).\cr
\code{dev_resids} returns a vector of squared (unscaled) deviance residuals (the summands defined for GLMs in McCullagh and Nelder 1989, p. 34; see Details of \code{\link{LL-family}} for other response families; these are also the squares of residuals returned by \code{\link[spaMM]{residuals}(., "deviance")}.). \cr
\code{deviance} returns the sum of squares of these (unscaled) deviance residuals, that is (consistently with \code{stats::deviance}) the unscaled deviance.\cr
\code{fitted} extracts fitted values.\cr
\code{response} extracts the response (as a vector).\cr
\code{fixef} extracts the fixed effects coefficients, \eqn{\beta}.\cr
\code{ranef} extracts the predicted random effects, \bold{Lv} (default since version 1.12.0), or optionally \bold{u} (see \code{\link{random-effects}} for definitions). \code{print.ranef} controls their printing. \cr
\code{getDistMat} returns a distance matrix for a geostatistical (Matérn etc.) random effect.\cr
\code{df.residual} extracts residual degrees-of-freedom for fitted models (here number of observations minus number of parameters of the model except residual dispersion parameters).
\code{wweights} extracts prior weights (as defined by the fitting functions's \code{prior.weights} argument).
}
\usage{
\method{formula}{HLfit}(x, which="hyper", ...)
\method{family}{HLfit}(object, ...)
\method{terms}{HLfit}(x, ...)
\method{nobs}{HLfit}(object, ...)
\method{logLik}{HLfit}(object, which, ...)
\method{fitted}{HLfit}(object, ...)
\method{fixef}{HLfit}(object, ...)
\method{ranef}{HLfit}(object, type = "correlated", ...)
\method{print}{ranef}(x, max.print = 40L, ...)
\method{deviance}{HLfit}(object, ...)
\method{df.residual}{HLfit}(object, ...)
\method{weights}{HLfit}(object, type, ...)
##
getDistMat(object, scaled=FALSE, which = 1L)
response(object,...)
dev_resids(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of class \code{HLfit}, as returned by the fitting functions in \code{spaMM}.}
  \item{type}{For \code{ranef}, use \code{type="correlated"} (default) to display the correlated random effects (\bold{Lv}), whether in a spatial model, or a random- coefficient model. Use \code{type="uncorrelated"} to pretty-print the elements of the \code{<object>$ranef} vector (\bold{u}).\cr
  For \code{weights}, either \code{"prior"} or \code{"working"}, with the same meaning as for \code{\link{weights.glm}}: respectively the prior weights, or the weights used in the final iteration of the IRLS algorithm.
  }
  \item{which}{* For \code{logLik}, the name of the element of the \code{APHLs} list to return (see Details for any further possibility). The default depends on the fitting method. In particular, if it was REML or one of its variants, the default is to return the log restricted likelihood (exact or approximated).\cr
  * For \code{getDistMat}, an integer, to select a random effect from several for which a distance matrix may be constructed.\cr
  * For \code{formula}, by default the model formula with non-expanded \code{multIMRF} random-effect terms is returned, while for \code{which=""} a formula with \code{multIMRF} terms expanded as \code{IMRF} terms is returned.}
  \item{scaled}{If \code{FALSE}, the function ignores the scale parameter \eqn{rho} and returns unscaled distance.}
  \item{x}{For \code{print.ranef}: the return value of \code{ranef.HLfit}.}
  \item{max.print}{Controls \code{options("max.print")} locally.}
  \item{\dots}{Other arguments that may be needed by some method.}
}
\value{
\code{formula} returns a \code{formula}, except a list of them from \code{fitmv()} output.

\code{terms} returns an object of class \code{c("terms", "formula")} which contains the \emph{terms} representation of a symbolic model.  See \code{\link[stats]{terms.object}} for its structure. \code{terms(<fitmv() result>)} returns a list of such terms.

Other return values are numeric (for \code{logLik}), vectors (most cases), matrices or dist objects (for \code{getDistMat}), or a family object (for \code{family}). \code{ranef} returns a list of vectors or matrices (the latter for random-coefficient terms). 
}
\seealso{
See \code{\link{summary.HLfit}} whose return value include the tables of fixed-effects coefficients and random-effect variances displayed by the summary, \code{\link{residuals.HLfit}} to extract various residuals, \code{\link{residVar}} to extract residual variances or information about residual variance models, \code{\link{hatvalues}} to extract leverages, \code{\link{get_matrix}} to extract the model matrix and derived matrices, and \code{\link{vcov.HLfit}} to extract covariances matrices from a fit, \code{\link{get_RLRsim_args}} to extract arguments for (notably) tests of random effects in LMMs.
}
\references{
McCullagh, P. and Nelder J. A. (1989) Generalized linear models. Second ed. Chapman & Hall: London.
}
\examples{
data("wafers")
m1 <- fitme(y ~ X1+X2+(1|batch), data=wafers)
fixef(m1)
ranef(m1)

data("blackcap")
fitobject <- fitme(migStatus ~ 1 + Matern(1|longitude+latitude),data=blackcap,
                       fixed=list(nu=4,rho=0.4,phi=0.05))
getDistMat(fitobject)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
