\name{options}
\alias{spaMM.options}
\alias{spaMM.getOption}
%~ \alias{.SpaMM}
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
  \item{x}{a character string holding an option name.}
  \item{\dots}{A named value or a list of named values. The following values, with their defaults, 
    are used in \code{spaMM}:

  \describe{
\item{\code{ff_threshold}:}{Minimal number of matrix elements for some matrices to be stored on disk (using the \code{ff} package) rather than in memory}
  

\item{\code{COMP_maxn}:}{Number of terms for truntion of infinite sums that are evaluated in the fitting of \code{\link{COMPoisson}} models.}

%\item{\code{RHOMAX}:}{Maximum value of the scale parameter \eqn{\rho} of the Matérn correlation function. Reserved variable, not set by default, but set by \code{corrHLfit}.}

%\item{\code{NUMAX=50}:}{Maximum value of the smoothness parameter  \eqn{\nu} of the Matérn correlation function}

\item{\code{TRACE.UNLINK=FALSE}:}{Whether to delete preexisting files named as \code{HLCor.args.*.RData} when the trace option is used in corrHLfit or HLCor.}

\item{\code{MESSAGES.FULL.STACK=TRUE}:}{Whether to give information on the stack of calls in some warning messages.}

\item{\code{QRmethod}:}{To control the methods used in intensive matrix computations. Computations may be faster in simple block random effect models when this is set to \code{"Matrix::qr"}.}

\item{\code{USEEIGEN=TRUE}:}{Whether to use the Eigen C++ library for some matrix computations.}

\item{\code{USElmwithQ=FALSE}:}{Whether to use the Eigen C++ library for solving linear systems.}

\item{\code{USEEIGEN_QR=TRUE}:}{Whether to use the Eigen C++ library for some other QR computations.}

\item{\code{maxLambda=1e10}:}{The maximum value of lambda: higher fitted 
 lambda values in HLfit are reduced to this.} 
\item{\code{example_maxtime=2.5}:}{Used in the documentation to control whether the longer examples should be run. 
The approximate running time of given examples on one author's laptop is compared to this value.} 
}
and possibly other undocumented values for development purposes.

} % \item{\dots}{}
} % arguments
\details{
   
  Invoking \code{spaMM.options()} with no arguments returns a list with the
  current values of the options. Invoking \code{spaMM.getOption(<option name>)} returns the value of the option rather than a list.

%~  Internally, the options are kept in the list \code{.SpaMM}.\cr

  spaMM.options() provides an interface for changing maximal values of parameters of the Matérn correlation function. However, 
  it is not recommended to change these values unless a spaMM message specifically suggests so. Errors may occur if too low values are chosen as these may conflit with default initial values for the parameters.
}
\value{
  For \code{spaMM.getOption}, the current value set for option \code{x}, or
  \code{NULL} if the option is unset.

  For \code{spaMM.options()}, a list of all set options sorted by category.  For
  \code{spaMM.options(name)}, a list of length one containing the set value,
  or \code{NULL} if it is unset.  For uses setting one or more options,
  a list with the previous values of the options changed (returned
  invisibly).
}
\examples{
  spaMM.options()
  spaMM.getOption("example_maxtime")
  \dontrun{
  spaMM.options(maxLambda=1e06)
  }
}

