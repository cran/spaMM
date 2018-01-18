\name{summary.HLfit}
\alias{summary}
\alias{summary.HLfit}
\alias{summary.HLfitlist}
\alias{summary.fixedLRT}
\alias{print}
\alias{print.HLfit}
\alias{print.HLfitlist}
\alias{print.fixedLRT}

\title{
  Summary and print methods for fit and test results.
}
\description{
  Summary and print methods for results from HLfit or related functions.  
}
\usage{
\method{summary}{HLfit}(object, details=FALSE, ...)
\method{summary}{HLfitlist}(object, ...)
\method{summary}{fixedLRT}(object, verbose=TRUE, ...)
\method{print}{HLfit}(x,...)
\method{print}{HLfitlist}(x,...)
\method{print}{fixedLRT}(x,...)
}
\arguments{
  \item{object}{
  The return object of HLfit or related functions.
}
  \item{x}{
  The return object of HLfit or related functions.
}
  \item{verbose}{
  for \code{summary.fixedLRT}, whether to print the model fits or not.
}
  \item{details}{
  whether to print some obscure details. Currently affects only random-coefficients models (see Details).
}
  \item{\dots}{
  further arguments passed to or from other methods. 
}
}
\value{
  These methods return the object invisibly. They print details of the (lower level) HLfit results in a convenient form.
}
\details{The random effect terms of the linear predictor are of the form \bold{ZLv}. In particular, for \bold{random-coefficients models} (i.e., including random-effect terms such as \code{(z|group)} specifying a random-slope component), correlated random effects are represented as \bold{b = Lv} for some matrix \bold{L}, and where the elements of \bold{v} are uncorrelated. In the output of the fit, the \code{Var.} column gives the
variances of the correlated effects, \bold{b=Lv}. The \code{Corr.} column(s) give their correlation(s). If \code{details} is TRUE, estimates and SEs of the (log) variances of the elements of \bold{v} are reported as for other random effects in the \code{Estimate} and \code{cond.SE.} columns of the table of lambda coefficients. However, this non-default output is potentially misleading as the elements of \bold{v} cannot generally be assigned to specific terms (such as intercept and slope) of the random-effect formula, and the representation of \bold{b} as \bold{Lv} is not unique.
}
\examples{
## see examples of corrHLfit usage
}
\keyword{ print }
