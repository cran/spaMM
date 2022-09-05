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
  Summary and print methods for results from HLfit or related functions. \code{summary} may also be used as an extractor (see e.g. \code{\link{beta_table}}).  
}
\usage{
\method{summary}{HLfit}(object, details=FALSE, max.print=100L, verbose=TRUE, ...)
\method{summary}{HLfitlist}(object, ...)
\method{summary}{fixedLRT}(object, verbose=TRUE, ...)
\method{print}{HLfit}(x,...)
\method{print}{HLfitlist}(x,...)
\method{print}{fixedLRT}(x,...)
}
\arguments{
  \item{object}{An object of class \code{HLfit}, as returned by the fitting functions in \code{spaMM}.}
  \item{x}{  The return object of HLfit or related functions.}
  \item{verbose}{ For \code{summary.HLfit}, whether to print the screen output that is the primary purpose of summary. \code{verbose=FALSE} may be convenient when \code{summary} is used as an extractor. For \code{summary.fixedLRT}, whether to print the model fits or not. }
  \item{max.print}{Controls \code{options("max.print")} locally.}
  \item{details}{A vector with elements controlling whether to print some obscure details. Element \code{ranCoefs=TRUE} will print details about random-coefficients terms (see Details); and element \code{p_value="Wald"} will print a p-value for the t-value of each fixed-effect coefficient, assuming a gaussian distribution of the test statistic (but, beyond the generally questionable nature of p-value tables, see e.g. \code{LRT} and \code{fixedLRT} for alternative testing approaches).}
  \item{\dots}{  further arguments passed to or from other methods. }
}
\value{
  The return value is a list whose elements may be subject to changes, but two of them can be considered stable, and are thus part of the API: the \code{beta_table} and \code{lambda_table} which are the displayed tables for the coefficients of fixed effects and random-effect variances.
}
\details{The random effect terms of the linear predictor are of the form \bold{ZLv}. In particular, for \bold{random-coefficients models} (i.e., including random-effect terms such as \code{(z|group)} specifying a random-slope component), correlated random effects are represented as \bold{b = Lv} for some matrix \bold{L}, and where the elements of \bold{v} are uncorrelated. In the output of the fit, the \code{Var.} column gives the
variances of the correlated effects, \bold{b=Lv}. The \code{Corr.} column(s) give their correlation(s). If \code{details} is TRUE, estimates and SEs of the (log) variances of the elements of \bold{v} are reported as for other random effects in the \code{Estimate} and \code{cond.SE.} columns of the table of lambda coefficients. However, this non-default output is potentially misleading as the elements of \bold{v} cannot generally be assigned to specific terms (such as intercept and slope) of the random-effect formula, and the representation of \bold{b} as \bold{Lv} is not unique.
}
\examples{
## see examples of fitme() or corrHLfit() usage
}
\keyword{ print }
