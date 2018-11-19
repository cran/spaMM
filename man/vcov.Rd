\name{vcov}
\alias{vcov}
\alias{vcov.HLfit}
\alias{Corr}
\title{
Extract covariance or correlation matrices from a fitted model object}
\description{
\code{vcov} returns the variance-covariance matrix of the fixed-effects coefficients.
\code{Corr} returns a correlation matrix of random effects.
}
\usage{
\method{vcov}{HLfit}(object, ...)
Corr(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A object of class \code{HLfit}, as returned by the fitting functions in \code{spaMM}.}
  %\item{which}{not documented, for development purposes}
  \item{\dots}{Other arguments that may be needed by some method.}
}
\value{
\code{vcov} returns a matrix. \code{Corr} returns a list, for the different random effect terms. For each random-effect term, the returned element is a non-trivial unconditional correlation matrix of the vector \dQuote{\bold{v}} of random effects (\bold{v} as defined in see Details of \code{\link{HLfit}}) for this term, if there is any such matrix. Otherwise the returned element is a information message.
}
\examples{
data("wafers")
m1 <- HLfit(y ~X1+X2+(1|batch),
          resid.model = ~ 1 ,data=wafers,HLmethod="ML")
vcov(m1)
}
