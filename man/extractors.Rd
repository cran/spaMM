\name{extractors}
\alias{extractors}
\alias{predictionCoeffs}
\alias{getDistMat}
\alias{logLik}
\alias{logLik.HLfit}
\alias{fitted}
\alias{fitted.HLfit}
\alias{fixef}
\alias{fixef.HLfit}
\alias{ranef}
\alias{ranef.HLfit}
\alias{vcov}
\alias{vcov.HLfit}
\alias{Corr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Functions to extract various components of a fit
}
\description{
\code{logLik} extracts the log-likelihood (exact or approximated).
\code{fitted} extracts fitted values (see \code{\link{fitted.values}}).
\code{fixef} extracts the fixed effects coefficients, \eqn{\beta}.
\code{ranef} extracts the predicted random effects, \eqn{u}.
\code{vcov} returns the variance-covariance matrix of the fixed-effects coefficients.
\code{Corr} returns a correlation matrix of random effects (with restrictions, see Details).
\code{predictionCoeffs} extracts or precomputes coefficients for prediction (see \code{\link{predict}} for an example)
\code{getDistMat} extracts a distance matrix for a Matérn correlation model.
}
\usage{
\method{logLik}{HLfit}(object,REML,...)
\method{fitted}{HLfit}(object,...)
\method{fixef}{HLfit}(object,...)
% \method{ranef}{HLfit}(object,ranef.class=TRUE,...)
\method{ranef}{HLfit}(object,...)
\method{vcov}{HLfit}(object,...)
predictionCoeffs(object)
getDistMat(object,scaled=FALSE)
Corr(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{The return object of an HLfit or similar function.}
  \item{REML}{If \code{TRUE}, the function returns the log restricted likelihood (exact or approximated).}
  \item{scaled}{If \code{FALSE}, the function ignores the scale parameter \eqn{rho} and returns unscaled distance.}
%  \item{ranef.class}{If \code{TRUE}, the return value inherits from class \code{ranef}, so that the \code{print.ranef} method applies to it.}
  \item{\dots}{Other arguments that may be needed by some method.}
}
\value{
Return values are numeric (for \code{logLik}), vectors (most cases), matrices (for \code{vcov}), matrices or dist objects (for \code{getDistMat}). \code{ranef} returns a vector with attributes, 
which inherits from class \code{ranef} which has its own (undocumented) \code{print} method. 

\code{Corr} currently returns the correlation matrix of the random effects which are described as \bold{Lv} (see \code{\link{HLfit}})
}
\examples{
data(wafers)
m1 <- HLfit(y ~X1+X2+(1|batch),
          resid.formula = ~ 1 ,data=wafers,HLmethod="ML")
fixef(m1)
vcov(m1)
ranef(m1)
## see 'predict' for an example with predictionCoeffs and getDistMat
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
