\name{extractors}
\alias{extractors}
\alias{predictionCoeffs}
\alias{fitted}
\alias{fitted.HLfit}
\alias{fixef}
\alias{fixef.HLfit}
\alias{ranef}
\alias{ranef.HLfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Functions to extract various components of a fit
}
\description{
\code{fitted} extracts fitted values (see \code{\link{fitted.values}}).
\code{fixef} extracts the fixed effects coefficients, \eqn{\beta}.
\code{ranef} extracts the predicted random effects, \eqn{u}.
\code{predictionCoeffs} precomputes coefficients for prediction (see \code{\link{predict}} for an example)
}
\usage{
\method{fitted}{HLfit}(object,...)
\method{fixef}{HLfit}(object,...)
\method{ranef}{HLfit}(object,...)
predictionCoeffs(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{The return object of an HLfit or similar function.}
  \item{\dots}{Other arguments that may be needed by some method.}
}
\value{
\code{predictionCoeffs} returns a numeric vector.
}
\examples{
## see 'predict' for a example with predictionCoeffs 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
