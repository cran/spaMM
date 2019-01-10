\name{corr_family}
\alias{corr_family}
\alias{SAR_WWt}
\alias{print.corr_family}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
\code{corr_family} objects
}
\description{
\code{corr_family} objects provide a convenient way to implement correlation models handled by \code{spaMM}, analogous to \code{family} objects. These objects are undocumented (but there are documentation pages for each of the models implemented).
}
\usage{
# Matern(...)           # see help(Matern)
# Cauchy(...)           # see help(Cauchy)
# corrMatrix(...)       # see help(corrMatrix)
# AR1(...)              # see help(AR1)
# adjacency(...)        # see help(adjacency)
% MRF(...)              # see help(MRF)    % F I X M E
\method{print}{corr_family}(x,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{\code{corr_family} object.}
  \item{\dots}{arguments that may be needed by some \code{corr_family} object or some print method.}
}
\keyword{ models }% use one of  RShowDoc("KEYWORDS")
\keyword{ spatial }% __ONLY ONE__ keyword per line
