\name{options}
\alias{spaMM.options}
\alias{spaMM.getOption}
\alias{.SpaMM}
\title{spaMM options settings}
\description{
  Allow the user to set and examine a variety of \emph{options}
  which affect operations of the spaMM package.
}
\usage{
spaMM.options(\dots)

spaMM.getOption(x)

}
\arguments{
  \item{\dots}{Anything consistent with stadard format for function arguments. However, only the following ones
    are used in \code{spaMM}:

  \describe{

\item{\code{RHOMAX=100000}:}{Maximum value of the scale parameter  \eqn{\rho} of the Matern correlation function}

\item{\code{NUMAX=50}:}{Maximum value of the smoothness parameter  \eqn{\nu} of the Matern correlation function}

\item{\code{TRACE.UNLINK=FALSE}:}{Whether to delete preexisting files named as \code{HLCor.args.*.RData} when the 
       trace option is used in corrHLfit or HLCor.}

\item{\code{INIT.HLFITNAME=NA}:}{Something users should not care about.}

\item{\code{TMPNAME=NA}:}{Something users should not care about.}

}
  }

  \item{x}{a character string holding an option name.}
}
\details{
   
  Invoking \code{spaMM.options()} with no arguments returns a list with the
  current values of the options. Invoking \code{spaMM.getOption(<option name>)} returns the value of the option rather than a list.

  Internally, the options are kept in the list \code{.SpaMM}.\cr

  spaMM.options() provides an interface for changing maximal values of parameters of the Matern correlation function. However, 
  currently it is not recommended to change these values. Errors may occur if too low valeus are chosen.
}
\value{
  For \code{spam.getOption}, the current value set for option \code{x}, or
  \code{NULL} if the option is unset.

  For \code{spaMM.options()}, a list of all set options sorted by category.  For
  \code{spaMM.options(name)}, a list of length one containing the set value,
  or \code{NULL} if it is unset.  For uses setting one or more options,
  a list with the previous values of the options changed (returned
  invisibly).
}
\examples{
  spaMM.options()
}

