\name{summary.HLfit}
\alias{summary}
\alias{summary.HLfit}
\alias{summary.HLCor}
\alias{summary.corrHLfit}
\alias{summary.fixedLRT}
\alias{print}
\alias{print.HLfit}
\alias{print.HLCor}
\alias{print.corrHLfit}
\alias{print.fixedLRT}

\title{
  Summary and print methods for fit and test results.
}
\description{
  Summary and print methods for results from HLfit / HLCor / corrHLfit / fixedLRT.  
}
\usage{
\method{summary}{HLfit}(object,...)
\method{summary}{HLCor}(object,...)
\method{summary}{corrHLfit}(object,...)
\method{summary}{fixedLRT}(object,verbose=TRUE,...)
\method{print}{HLfit}(x,...)
\method{print}{HLCor}(x,...)
\method{print}{corrHLfit}(x,...)
\method{print}{fixedLRT}(x,...)
}
\arguments{
  \item{object}{
  The return object of an HLCor / HLfit / corrHLfit call.
}
  \item{x}{
  The return object of an HLCor / HLfit / corrHLfit call.
}
  \item{verbose}{
  for \code{summary.fixedLRT}, whether to print the model fits or not.
}
  \item{\dots}{
  further arguments passed to or from other methods. 
}
}
\value{
  These methods return the object invisibly. They print details of the (lower level) HLfit results in a convenient form.
}
\examples{
## see examples of corrHLfit usage
}
\keyword{ print }
