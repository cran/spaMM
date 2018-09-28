\name{how}
\alias{how}
\alias{how.default}
\alias{how.HLfit}
\title{
Extract information about how an object was obtained
}
\description{
\code{how} is defined as a generic with currently only one non-default method, for objects of class \code{HLfit}. This method provide information about how such a fit was obtained.
}
\usage{
how(object, ...)
\method{how}{HLfit}(object, devel=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{Any R object.}
  \item{devel}{Boolean; Whether to provide additional cryptic information. For development purposes, not further documented.}
  \item{\dots}{Other arguments that may be needed by some method.}
}
\value{
  A list, returned invisibly, whose elements are not further described here. The function prints a message presenting these elements, some of which may be slightly cryptic. This function is work in progress.
}
\examples{
foo <- HLfit(y~x, data=data.frame(x=runif(3),y=runif(3)),HLmethod="ML",ranFix=list(phi=1))
how(foo)
}
