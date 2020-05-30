\name{options}
\alias{spaMM.options}
\alias{spaMM.getOption}
\alias{LevenbergM}
\alias{barstyle}
\title{spaMM options settings}
\description{
  Allow the user to set and examine a variety of \emph{options}
  which affect operations of the spaMM package.
}
\usage{
spaMM.options(\dots, warn = TRUE)

spaMM.getOption(x)

}
\arguments{
  \item{x}{a character string holding an option name.}
  \item{warn}{Boolean: whether to warn if a previously undefined options is being defined (a protection against typos).}
  \item{\dots}{A named value or a list of named values. The following values, with their defaults, 
    are used in \code{spaMM}:

  \describe{
\item{\code{LevenbergM=NULL}:}{NULL or boolean. Whether to use a Levenberg-Marquardt-like algorithm (see Details) by default in most computations. But it is advised to use instead \code{control.HLfit=list(LevenbergM=...)} to control this on a case-by-case basis. The joint default behaviour is that Levenberg-Marquardt is used by default for binomial response data that takes only extreme values (in particular, for binary 0/1 response), and that for other models the fitting algorithm switches to it if divergence is suspected. \code{FALSE} inhibits its use; \code{TRUE} forces its use for all iterative least-square fits, except when 'confint()' is called.}

\item{\code{example_maxtime=0.7}:}{Used in the documentation to control whether the longer examples should be run. 
The approximate running time of given examples on one author's laptop is compared to this value.} 
\item{\code{optimizer1D="optimize"}:}{Optimizer for one-dimensional optimization. If you want to control the initial value, you should select another optimizer.}
\item{\code{optimizer=".safe_opt"}:}{Optimizer for optimization in several dimensions. Use \code{optimizer="nloptr"} to call \code{\link[nloptr]{nloptr}} with method \code{"NLOPT_LN_BOBYQA"}; use \code{optimizer="bobyqa"} to call \code{\link[minqa]{bobyqa}}; and use \code{optimizer="L-BFGS-B"} to call \code{\link{optim}} with method \code{"L-BFGS-B"}. The default \code{".safe_opt"} uses \code{nloptr} except in some cases where it expects or detects problems with it (the source coe should be consulted for details). The optimizer can also be specified on a fit-by-fit basis as the value of \code{control$optimizer} in a \code{fitme} call, or as the value of \code{control.corrHLfit$optimizer}.
}
\item{\code{nloptr}:}{Default control values of \code{nloptr} calls.} 
\item{\code{bobyqa}:}{Default control values of \code{bobyqa} calls.} 

\item{\code{maxLambda=1e10}:}{The maximum value of lambda: higher fitted lambda values in HLfit are reduced to this. Since version 3.1.0, a much smaller lambda bound is deduced from \code{maxLambda} for COMPoisson and log-link response families.}
\item{\code{regul_lev_lambda}}{Numeric (default: 1e-8); lambda leverages numerically 1 are replaced by 1- \code{regul_lev_lambda}}  
\item{\code{COMP_maxn}:}{Number of terms for truncation of infinite sums that are evaluated in the fitting of \code{\link{COMPoisson}} models.}
\item{\code{CMP_asympto_cond}:}{Condition for applying an approximation or the COMPoisson response family, as detailed in \code{\link{COMPoisson}}.} 
\item{\code{Gamma_min_y=1e-10}:}{A minimum response value in Gamma-response models; used to check data, and in \code{simulate()} to correct the simulation results.} 

\item{\code{QRmethod}:}{A character string, to control whether dense matrix or sparse matrix methods are used in intensive matrix computations, overcoming the defaults choices made by \code{spaMM} in this respect. Possible values are \code{"dense"} and \code{"sparse"}.}

\item{\code{matrix_method}:}{A character string, to control the factorization of dense model matrices. Default value is \code{"def_sXaug_EigenDense_QRP_scaled"}. The source code should be consulted for further information.}

\item{\code{Matrix_method}:}{A character string, to control the factorization of sparse model matrices. Default value is \code{"def_sXaug_Matrix_QRP_scaled"}. The source code should be consulted for further information.}

\item{\code{stylefns}:}{Default colors of some screen output (notably that of some fitting functions when called with argument \code{verbose=c(TRACE=TRUE)})}
\item{\code{barstyle}:}{Integer, or Boolean interpreted as Integer, or quoted expression evaluating to such types; controlling the display of some progress bars. If zero, no progress bar should be displayed; otherwise, a bar should be displayed. Further, when \code{\link[utils]{txtProgressBar}} is called, \code{barstyle} is passed as its \code{style} argument. Default is \code{quote(if(interactive()) {3L} else {0L})}.}
}
and many other undocumented values for programming or development purposes. Additional options without default values can also be used (e.g., see \code{\link{sparse_precision}}).

} % \item{\dots}{}
} % arguments
\details{
   
  spaMM.options() provides an interface for changing maximal values of parameters of the Matérn correlation function. However, 
  it is not recommended to change these values unless a spaMM message specifically suggests so.
  
  By default spaMM use Iteratively Reweighted Least Squares (IRLS) methods to estimate fixed-effect parameters (jointly with predictions of random effects). However, a Levenberg-Marquardt algorithm, as described by Nocedal & Wright (1999, p. 266), is also implemented. The Levenberg-Marquardt algorithm is designed to optimize a single objective function with respect to all its parameters. It is thus well suited to compute a PQL fit, which is based on maximization of a single function, the h-likelihood. By contrast, in a fit of a mixed model by (RE)ML, one computes jointly fixed-effect estimates that maximizes marginal likelihood, and random-effect values that maximize h-likelihood given the fixed-effect estimates. The gradient of marginal likelihood with respect to fixed-effect coefficients does not generally vanishes at the solution (although it remains close to zero except in \dQuote{difficult} cases with typically little information in the data). The Levenberg-Marquardt algorithm is not directly applicable in this case, as it may produce random-effect values that increases marginal likelihood rather than h-likelihood. The (RE)ML variant of the algorithm implemented in spaMM may therefore use additional nested h-likelihood-maximizing steps for correcting random-effect values. In version 3.1.0 this variant was revised for improved performance in difficult cases. 

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

