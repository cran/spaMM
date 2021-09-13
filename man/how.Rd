\name{how}
\alias{how}
\alias{how.default}
\alias{how.HLfit}
\alias{how.HLfitlist}
\title{
Extract information about how an object was obtained
}
\description{
\code{how} is defined as a generic with currently only one non-default method, for objects of class \code{HLfit}. This method provide information about how such a fit was obtained.
}
\usage{
how(object, ...)
\method{how}{HLfit}(object, devel=FALSE, verbose=TRUE, format=print, ...)
\method{how}{HLfitlist}(object, devel=FALSE, verbose=TRUE, format=print, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{Any R object.}
  \item{devel}{Boolean; Whether to provide additional cryptic information. For development purposes, not further documented.}
  \item{verbose}{Boolean; Whether to print information about the input object.}
  \item{format}{wrapper for printing format. E.g., \code{cat(crayon::yellow(s),"\n")} could be used instead of the default.}
  \item{\dots}{Other arguments that may be needed by some method.}
}
\value{
  A list, returned invisibly, whose elements are not further described here, some being slightly cryptic or subject to future changes However, \code{how(.)$fit_time} is a clean way of getting the fit time. If \code{verbose} is \code{TRUE}, the function prints a message presenting some of these elements.
}
\examples{
foo <- HLfit(y~x, data=data.frame(x=runif(3), y=runif(3)), method="ML", ranFix=list(phi=1))
how(foo)
}
