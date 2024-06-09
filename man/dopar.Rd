\name{dopar}
\alias{dopar}
\alias{combinepar}
\title{
Interface for parallel computations
}
\description{
\code{dopar} and \code{combinepar} are interfaces primarily designed to apply some function \code{fn} in parallel on columns of a matrix, although other uses are possible. Depending on the \code{nb_cores} argument, parallel or serial computation is performed. A socket cluster is used by default for parallel computations, but a fork cluster can be requested on linux and alike operating systems by using argument \code{cluster_args=list(type="FORK")}.

\code{dopar} has been designed to provide by default a progress bar in all evaluations contexts. A drawback is that different procedures are called depending e.g. on the type of cluster, with different possible controls. In particular, \code{foreach} is called in some cases but not others, so non-trivial values of its \code{.combine} control are not always enforced. The alternative interface \code{combinepar} will always use \code{foreach}, and will still try to provide by default a progress bar but may fail to do so in some cases (see Details).     

}
\usage{
dopar(newresp, fn, nb_cores = NULL, fit_env, 
      control = list(.final=function(v) if( ! is.list(v[[1]])) {do.call(cbind,v)} else v), 
      cluster_args = NULL, debug. = FALSE, iseed = NULL, 
      showpbar = eval(spaMM.getOption("barstyle")), 
      pretest_cores =NULL, ...)
combinepar(newresp, fn, nb_cores = NULL, cluster=NULL, fit_env, 
      control = list(.final=function(v) if( ! is.list(v[[1]])) {do.call(cbind,v)} else v), 
      cluster_args = NULL, debug. = FALSE, iseed = NULL, 
      showpbar = eval(spaMM.getOption("barstyle")), 
      pretest_cores =NULL, ...)

}
\arguments{
  \item{newresp}{
A matrix on whose columns \code{fn} will be applied (e.g., as used internally in \pkg{spaMM}, the return value of a \code{simulate.HLfit()} call); or an integer, then converted to a trivial matrix \code{matrix(seq(newresp),ncol=newresp,nrow=1)}.
}
  \item{fn}{
Function whose first argument is named \code{y}. The function will be applied for \code{y} taken to be each column of \code{newresp}.
}
  \item{nb_cores}{
Integer. Number of cores to use for parallel computations. If >1 (and no cluster is provided by the \code{cluster} argument), a cluster of \code{nb_cores} nodes is created, used, and stopped on completion of the computation. Otherwise, no parallel computation is performed.
}
  \item{cluster}{
(for \code{combinepar} only): a cluster object (as returned by \code{parallel::makeCluster} or \code{parallel::makeForkCluster}). If this is used, the \code{nb_cores} and \code{cluster_args} arguments are ignored. The cluster is not stopped  on completion of the computation
}
  \item{fit_env}{
(for socket clusters only:)  An environment, or a list, containing variables to be exported on the nodes of the cluster (by \code{parallel::clusterExport}); e.g., \code{list(bar=bar)} to pass object \code{bar} to each node. The argument \code{control(.errorhandling = "pass")}, below, is useful to find out missing variables.
}
  \item{control}{
 A list following the \code{foreach} control syntax, even if \code{foreach} is not used. There are limitations when \code{dopar} (but not \code{combinepar}) is used, in all but the first case below:
 \enumerate{
 \item for socket clusters, with \code{doSNOW} attached, \code{foreach} is called with default arguments including
 \code{i = 1:ncol(newresp), .inorder = TRUE, .errorhandling = "remove", .packages = "spaMM"}, and further arguments taken from the present function's \code{control} argument, which may also be used to override the defaults. For example, \code{.errorhandling = "pass"} is useful to get error messages from the nodes, and therefore \bold{strongly recommended} when first experimenting with this function.     
 \item for socket clusters, with \code{doSNOW} \bold{not} attached, \code{dopar} calls \code{pbapply} instead of \code{foreach} but \code{control$.packages} is still handled. The  result is still in the format returned in the first case, i.e. by \code{foreach}, taking the \code{control} argument into account. \code{pbapply} arguments may be passed through the \dots argument. 
 \item if a fork cluster is used, \code{dopar} calls \code{mclapply} instead of \code{foreach}. \code{control$mc.silent} can be used to control the \code{mc.silent} argument of \code{\link[parallel]{mclapply}}.
 \item (if \code{nb_cores=1} \code{dopar} calls \code{mclapply}).
 }
}
  \item{cluster_args}{
 A list of arguments passed to \code{parallel::makeCluster}. E.g., \code{outfile="log.txt"} may be useful to collect output from the nodes, and \code{type="FORK"} to force a fork cluster on linux(-alikes).
}
  \item{debug.}{
(for socket clusters only:)   For debugging purposes. Effect, if any, is to be defined by the \code{fn} as provided by the user.
}
  \item{iseed}{
(all parallel contexts:) Integer, or NULL. If an integer, it is used as the \code{iseed} argument of \code{\link[parallel]{clusterSetRNGStream}} to initialize \code{"L'Ecuyer-CMRG"} random-number generator (see Details). If \code{iseed} is \code{NULL}, the default generator is selected on each node, where its seed is not controlled. 
}
  \item{showpbar}{
(for socket clusters only:) Controls display of progress bar. See \code{\link{barstyle}} option for details.
}
  \item{pretest_cores}{
(for socket clusters only:) A function to run on the cores before running \code{fn}. It may be used to check that all arguments of the \code{fn} can be evaluated in the cores' environments (the internal function \code{.pretest_fn_on_cores} provides an example).   
}
  \item{\dots}{
Further arguments to be passed (unevaluated) to \code{fn}, if not caught on the way by \code{pbapply} (which means that different results may in principle be obtained depending on the mode of parallelisation, which is the kind of design issues that \code{combinepar} aims to resolve by always calling \code{foreach}).
}
}
\value{
The result of calling \code{foreach}, \code{pbapply} or \code{mclapply}, as dependent on the \code{control} argument and the interface used. A side-effect of either interface is to show a progress bar whose character informs about the type of parallelisation performed: a \code{"F"} or default \code{"="} character for fork clusters, a \code{"P"} for parallel computation via \code{foreach} and \code{doSNOW}, a \code{"p"} for parallel computation via \code{foreach} and \code{doFuture} or via \code{pbapply}, and \code{"s"} for serial computation \code{foreach} and \code{doParallel} or via \code{pbapply}.
}
\details{
Control of random numbers through the \code{"L'Ecuyer-CMRG"} generator and the \code{iseed} argument is not sufficient for consistent results when the \code{doSNOW} parallel backend is used, so if you really need such control in a \code{fn} using random numbers, do not use \code{doSNOW}. Yet, it is fine to use \code{doSNOW} for bootstrap procedures in spaMM, because the fitting functions do not use random numbers: only sample simulation uses them, and it is not performed in parallel.
%  (alternative function \code{\link{dofuture}} has no doSNOW related code: why?) _F I X M E__

\code{combinepar} calls \code{foreach::\%dopar\%} which assumes that a cluster has been declared using a suitable backend such as \code{doSNOW}, \code{doFuture} or \code{doParallel}. If only the latter is available, no progress bar is displayed. A method to render a bar when \code{doParallel} is used can be found on the Web, but that bar is not a valid progress bar as it is displayed only after all the processes have been run.     
}
\seealso{
\code{\link{dofuture}} is yet another interface with (essentially) the same functionalities as \code{dopar}. See the documentation of the \code{\link{wrap_parallel}} option for its differences from \code{dopar}.
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

% cf also help(Infusion::SLRT)
}

