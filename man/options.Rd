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
  

\item{\code{COMP_maxn}:}{Number of terms for truncation of infinite sums that are evaluated in the fitting of \code{\link{COMPoisson}} models.}

\item{\code{MESSAGES.FULL.STACK=TRUE}:}{Whether to give information on the stack of calls in some warning messages.}

\item{\code{QRmethod}:}{A character string, to control whether dense matrix or sparse matrix methods are used in intensive matrix computations, overcoming the defaults choices made by \code{spaMM} in this respect. Possible values are \code{"dense"} and \code{"sparse"}.}

\item{\code{matrix_method}:}{A character string, to control the factorization of dense model matrices. Default value is \code{"def_sXaug_EigenDense_QR_scaled"}. \code{"def_sXaug_EigenDense_QRP_scaled"} is sometimes useful (see \code{\link{arabidopsis}}). The source code should be consulted for further information.}

\item{\code{Matrix_method}:}{A character string, to control the factorization of sparse model matrices. Default value is \code{"def_sXaug_Matrix_QRP_scaled"}. The source code should be consulted for further information.}

\item{\code{LevenbergM=FALSE}:}{Boolean. Whether to use a Levenberg-Marquardt algorithm (see Details) by default in most computations. But it is advised to use \code{control.HLfit=list(LevenbergM=...)} to control this on a case-by-case basis.}

\item{\code{USEEIGEN=TRUE}:}{Whether to use the Eigen C++ library for some matrix computations.. The source code should be consulted for further information.}

\item{\code{wRegularization=FALSE}:}{Whether to to warn about the use of regularization in some operations on nearly singular matrices.}

\item{\code{maxLambda=1e10}:}{The maximum value of lambda: higher fitted 
 lambda values in HLfit are reduced to this.} 
\item{\code{example_maxtime=0.8}:}{Used in the documentation to control whether the longer examples should be run. 
The approximate running time of given examples on one author's laptop is compared to this value.} 
}
and possibly other undocumented values for development purposes.

} % \item{\dots}{}
} % arguments
\details{
   
  spaMM.options() provides an interface for changing maximal values of parameters of the Matérn correlation function. However, 
  it is not recommended to change these values unless a spaMM message specifically suggests so. Errors may occur if too low values are chosen as these may conflit with default initial values for the parameters.
  
  The implemented Levenberg-Marquardt algorithm is the variant described by Nocedal & Wright (1999, p. 266). 
}
\references{
Jorge Nocedal and Stephen J. Wright (1999) Numerical Optimization. Springer-Verlag, New York.
}
\value{
  For \code{spaMM.getOption}, the current value set for option \code{x}, or
  \code{NULL} if the option is unset.

  For \code{spaMM.options()}, a list of all set options.  For
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

