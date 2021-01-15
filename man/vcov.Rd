\name{vcov}
\alias{beta_table}
\alias{vcov}
\alias{vcov.HLfit}
\alias{Corr}
\alias{VarCorr}
\alias{VarCorr.HLfit}
\title{
Extract covariance or correlation components from a fitted model object}
\description{
\code{summary(<fit object>)$beta_table} returns the table of fixed-effect coefficients as it is printed by \code{summary}, including standard errors and t-values.
\code{vcov} returns the variance-covariance matrix of the fixed-effects coefficients.
\code{Corr} returns a correlation matrix of random effects.

\code{VarCorr} returns (co)variance parameters of random effects, and optionally residual variance(s), from a fit object, in a data frame format roughly consistent with the method of objects of class \code{"lme"}, in particular including columns with consistent names for easier extraction. One may have to consult the summary of the object to check the meaning of the contents of this data frame (e.g., of 'variance' coefficients for non-gaussian random effects). Other extractors to consider are \code{\link{get_ranPars}} and \code{\link{get_inits_from_fit}}, the latter providing parameters in a form suitable for initializing a fit.
}
\usage{
\method{vcov}{HLfit}(object, ...)
\method{VarCorr}{HLfit}(x, sigma = 1, add_residVars=TRUE, ...)
Corr(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object, x}{A fitted model object, inheriting from class \code{"HLfit"}, as returned by the fitting functions in \code{spaMM}.}
  %\item{which}{not documented, for development purposes}
  \item{add_residVars}{
  Boolean; whether to include residual variance information in the returned table.
  }
    \item{sigma}{
ignored argument, included for consistency with the generic function. 
}
  \item{\dots}{Other arguments that may be needed by some method.}
}
\value{
\code{vcov} returns a matrix. \code{Corr} returns a list, for the different random effect terms. For each random-effect term, the returned element is a non-trivial unconditional correlation matrix of the vector \dQuote{\bold{v}} of random effects (\bold{v} as defined in see Details of \code{\link{HLfit}}) for this term, if there is any such matrix. Otherwise the returned element is a information message. 

\code{VarCorr} returns either NULL (if no variance to report, as for a poisson GLM) or a data frame with columns for the grouping factor, term, variance of random effect,  standard deviation (the root of the variance), and optionally for correlation of random effect in random-coefficient terms.  Information about the residual variance is optionally included as the last row(s) of the data frame, when relevant (gaussian- or Gamma-response models with single scalar parameter; beware the meaning of the residual variance parameter for Gamma-response models).
}
\seealso{
\code{\link{get_inits_from_fit}} and \code{\link{get_ranPars}}.
}
\examples{
data("wafers")
m1 <- HLfit(y ~ X1+X2+(1|batch), resid.model = ~ 1 ,data=wafers, method="ML")
vcov(m1)

# Example from VarCorr() documentation in 'nlme' package
data("Orthodont",package = "nlme")
sp1 <- fitme(distance ~ age+(age|Subject), data = Orthodont, method="REML")
VarCorr(sp1)
}
