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

\code{VarCorr} returns (co)variance parameters of random effects, and residual variance, from a fit object of class \code{"HLfit"}, in a format roughly consistent with the method of objects of class \code{"lme"}, in particular including columns with consistent names for easier extraction. However, the structure of the return value is not fully specified and some model structures (e.g., including \code{adjacency} terms) may not be correctly handled.
}
\usage{
\method{vcov}{HLfit}(object, ...)
\method{VarCorr}{HLfit}(x, sigma = 1, message. = TRUE, ...)
Corr(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object, x}{A fitted model object, inheriting from class \code{"HLfit"}, as returned by the fitting functions in \code{spaMM}.}
  %\item{which}{not documented, for development purposes}
    \item{sigma}{
ignored argument, included for consistency with the generic function. 
}
  \item{message.}{
boolean; whether to print a reminder about the function.
}
  \item{\dots}{Other arguments that may be needed by some method.}
}
\value{
\code{vcov} returns a matrix. \code{Corr} returns a list, for the different random effect terms. For each random-effect term, the returned element is a non-trivial unconditional correlation matrix of the vector \dQuote{\bold{v}} of random effects (\bold{v} as defined in see Details of \code{\link{HLfit}}) for this term, if there is any such matrix. Otherwise the returned element is a information message. 

\code{VarCorr} returns a data frame with columns for the grouping factor, term, variance of random effect,  standard deviation (the root of the variance), and optionnally for correlation of random effect in random-coefficient terms.  Information about the residual variance is included as the last row in the data frame, when relevant (gaussian- or Gamma-response models with single scalar parameter; beware the meaing of the residual variance parameter for Gamma-response models).
}
\examples{
data("wafers")
m1 <- HLfit(y ~X1+X2+(1|batch),
          resid.model = ~ 1 ,data=wafers,HLmethod="ML")
vcov(m1)

# Example from VarCorr() documentation in 'nlme' package
data("Orthodont",package = "nlme")
sp1 <- fitme(distance ~ age+(age|Subject), data = Orthodont, method="REML")
VarCorr(sp1)
}
