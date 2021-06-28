\name{dopar}
\alias{dopar}
\title{
Interface for parallel computations
}
\description{
An interface to apply some function \code{fn} in parallel on columns of a matrix. \code{dopar} is used internally in spaMM, but is not logically restricted to mixed-effect applications, hence it can be used more widely. Depending on the \code{nb_cores} argument, parallel or serial computation is performed. A socket cluster is used by default for parallel computations, but a fork cluster can be requested on linux and alike operating systems by using argument \code{cluster_args=list(type="FORK")}.
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
(for socket clusters only:)  An environment, or a list, containing variables to be exported on the nodes of the cluster (by \code{parallel::clusterExport}).
}
  \item{control}{
 A list. In particular if the \code{doSNOW} package is attached, \code{foreach} is called with default arguments including
 \code{i = 1:ncol(newresp), .combine = "cbind", .inorder = TRUE, .errorhandling = "remove", .packages = "spaMM"}.
\code{control} may be used to provide non-default values of these arguments. For example, \code{.errorhandling = "pass"} is useful to get error messages from the nodes. If \code{doSNOW} is not attached, the result is still in the format returned by \code{foreach} with default \code{.combine="cbind"} or possible non-default \code{.combine="rbind"}. 

Alternatively, if a for cluster is used, \code{control$mc.silent} can be used to control the \code{mc.silent} argument of \code{\link[parallel]{mclapply}}.
}
  \item{cluster_args}{
 A list of arguments passed to \code{parallel::makeCluster}. E.g., \code{outfile="log.txt"} may be useful to collect output from the nodes, and \code{type="FORK"} to force a fork cluster on linux(-alikes).
}
  \item{debug.}{
(for socket clusters only:)   For debugging purposes. Effect, if any, is to be defined by the \code{fn} as provided by the user.
}
  \item{iseed}{
(all parallel contexts:) Integer, or NULL. If an integer, it is used to initialize \code{"L'Ecuyer-CMRG"} random-number generator (\code{iseed} argument of \code{\link[parallel]{clusterSetRNGStream}}), yet not ensuring identical control of random number generator accross different parallelisation options (see alternative function \code{\link{dofuture}} for such control). If \code{iseed} is \code{NULL}, the default generator is selected on each node, where its seed is not controlled. 
}
  \item{showpbar}{
(for socket clusters only:) Controls display of progress bar. See \code{\link{barstyle}} option for details.
}
  \item{pretest_cores}{
(for socket clusters only:) A function to run on the cores before running \code{fn}. It may be used to check that all arguments of the \code{fn} can be evaluated in the cores' environments (the internal function \code{.pretest_fn_on_cores} provides an example).   
}
  \item{\dots}{
Further arguments to be passed (unevaluated) to \code{fn}.
}
}
\value{
The result of calling \code{foreach}, \code{pbapply} or \code{mclapply}, as dependent on the \code{control} argument. A side-effect of \code{dopar} is to show a progress bar that informs about the type of parallelisation performed: a default \code{"="} character for fork clusters, and otherwise\code{"P"} for parallel computation via \code{foreach} and \code{doSNOW}, \code{"p"} for parallel computation via \code{pbapply}, and \code{"s"} for serial computation  via \code{pbapply}.
}
\seealso{
\code{\link{dofuture}} for an alternative implementation of (essentially) the same functionalities, and \code{\link{wrap_parallel}} for its differences from \code{dopar}.
}
\examples{
## See source code of spaMM_boot()

\dontrun{
# Useless function, but requiring some argument beyond the first
foo <- function(y, somearg, ...) {
  if ( is.null(somearg) || TRUE ) length(y)
}

# Whether FORK can be used depends on OS and whether Rstudio is used:
dopar(matrix(1,ncol=4,nrow=3), foo, fit_env=list(), somearg=NULL, 
  nb_cores=2, cluster_args=list(type="FORK"))
}
}

