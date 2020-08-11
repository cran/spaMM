\name{dopar}
\alias{dopar}
\title{
Interface for parallel computations
}
\description{
An interface to apply some function \code{fn} in parallel on columns of a matrix. Depending on the \code{nb_cores} argument, parallel or serial computation is performed, and depending whether the \code{doSNOW} package is attached, \code{foreach} or \code{pbapply} is called (\code{doSNOW} allows more efficient load balancing than \code{pbapply}). It wraps these calls, handling all cluster management. It is used internally in spaMM, but is not logically restricted to mixed-effect applications, hence it can be used more widely. 
}
\usage{
dopar(newresp, fn, nb_cores = NULL, fit_env, control = list(), 
      cluster_args = NULL, debug. = FALSE, iseed = NULL, 
      showpbar = eval(spaMM.getOption("barstyle")), 
      pretest_cores =NULL, ...)
}
\arguments{
  \item{newresp}{
A matrix on whose columns \code{fn} will be applied (e.g., as used internally in spaMM, the return value of a \code{simulate.HLfit()} call). 
}
  \item{fn}{
Function, whose first argument is named \code{y}, to be applied to each column of \code{newresp}.
}
  \item{nb_cores}{
Integer. Number of cores to use for parallel computations. If >1, a cluster of \code{nb_cores} nodes is used. Otherwise, no parallel computation is performed.
}
  \item{fit_env}{
  A environment, or a list, containing variables to be exported on the nodes of the cluster (by \code{parallel::clusterExport}).
}
  \item{control}{
 A list. In particular if the \code{doSNOW} package is attached, \code{foreach} is called with default arguments including
 \code{i = 1:ncol(newresp), .combine = "cbind", .inorder = TRUE, .errorhandling = "remove", .packages = "spaMM"}.
\code{control} may be used to provide non-default values of these arguments. For example, \code{.errorhandling = "pass"} is useful to get error messages from the nodes. If \code{doSNOW} is not attached, the result is still in the format returned by \code{foreach} with default \code{.combine="cbind"} or possible non-default \code{.combine="rbind"}. 
}
  \item{cluster_args}{
 A list of arguments passed to \code{parallel::makeCluster}. E.g., \code{outfile="log.txt"} may be useful to collect output from the nodes.
}
  \item{debug.}{
  For debugging purposes. Effect, if any, is to be defined by the \code{fn} as provided by the user.
}
  \item{iseed}{
Integer, or NULL. If an integer, \code{parallel::clusterSetRNGStream(cl = cl, iseed)} is called, setting and initializing \code{"L'Ecuyer-CMRG"} random-number generator on the nodes. This is useful if random numbers are used on the nodes (see \code{\link[parallel:RngStream]{clusterSetRNGStream}}). If \code{iseed} is \code{NULL}, the default generator is selected on each node, where its seed is not controlled. 
}
  \item{showpbar}{
Controls display of progress bar. See \code{\link{barstyle}} option for details.
}
  \item{pretest_cores}{
A function to run on the cores before running \code{fn}. It may be used to check that all arguments of the \code{fn} can be evaluated in the cores' environments (the internal function \code{.pretest_fn_on_cores} provides an example).   
}
  \item{\dots}{
Further arguments to be passed (unevaluated) to \code{fn}.
}
}
\value{
The result of calling \code{foreach} or \code{pbapply}, as dependent on the \code{control} argument. A side-effect of \code{dopar} is to show a progress bar that informs about the type of parallelisation performed: \code{"P"} for parallel computation via \code{foreach} and \code{doSNOW}, \code{"p"} for parallel computation via \code{pbapply}, \code{"s"} for serial computation  via \code{pbapply}.
}
\examples{
## See source code of spaMM_boot()
}

