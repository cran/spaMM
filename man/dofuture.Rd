\name{dofuture}
\alias{dofuture}
\title{
Interface for parallel computations
}
\description{
interface to apply some function \code{fn} in parallel on columns of a matrix. It is not logically restricted to mixed-effect applications, hence it can be used more widely. Depending on the \code{nb_cores} argument, parallel or serial computation is performed, calling the \code{future.apply::future_apply} function. A socket cluster is used by default for parallel computations, but a fork cluster can be requested on linux and alike operating systems by using argument \code{cluster_args=list(type="FORK")}.
}
\usage{
dofuture(newresp, fn, nb_cores=NULL, fit_env, control=list(), 
      cluster_args=NULL, debug.=FALSE, iseed=NULL, 
      showpbar="ignored", pretest_cores="ignored",
      ... )
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
 A list. The only effective control is \code{.combine="rbind"} (mimicking the \code{foreach} syntax used in the alternative interface \code{\link{dopar}}). 
}
  \item{cluster_args}{
 A list of arguments passed to \code{parallel::makeCluster} or \code{parallel::makeForkCluster}. E.g., \code{outfile="log.txt"} may be useful to collect output from the nodes, and \code{type="FORK"} to force a fork cluster on linux(-alikes).
}
  \item{debug.}{
(for socket clusters only:)   For debugging purposes. Effect, if any, is to be defined by the \code{fn} as provided by the user.
}
  \item{iseed}{
Integer, or NULL. If an integer, it is used to initialize \code{"L'Ecuyer-CMRG"} random-number generator (\code{iseed} argument of \code{\link[parallel]{clusterSetRNGStream}}), with identical effect accross different models of parallelisation. If \code{iseed} is \code{NULL}, the seed is not controlled. 
}
  \item{showpbar,pretest_cores}{
Currently ignored; for consistency with \code{dopar} fromal arguments.
}
  \item{\dots}{
Further arguments to be passed (unevaluated) to \code{future.apply} (and then possibly to \code{fn}).
}
}
\value{
The result of calling \code{future.apply}.  If the \code{progressr} package is loaded, a side-effect of \code{dofuture} is to show a progress bar with character 'S' or 'P' or 'F' depending on parallelisation status (serial/socket/fork).
}
\seealso{
\code{\link{dopar}} for an alternative implementation of (essentially) the same functionalities, and \code{\link{wrap_parallel}} for its differences from \code{dofuture}.
}

\examples{

\dontrun{
if (requireNamespace("future.apply", quietly = TRUE)) {

 # Useless function, but requiring some argument beyond the first
 foo <- function(y, somearg, ...) {
   if ( is.null(somearg) || TRUE ) length(y)
 }

 # Whether FORK can be used depends on OS and whether Rstudio is used:
   dofuture(matrix(1,ncol=4,nrow=3), foo, fit_env=list(), somearg=NULL, 
     nb_cores=2, cluster_args=list(type="FORK"))
}
}
}
