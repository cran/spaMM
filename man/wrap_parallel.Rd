\name{wrap_parallel}
\alias{wrap_parallel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Selecting interfaces for parallelisation
}
\description{
spaMM implements two interfaces for parallelisation, \code{\link{dopar}} and \code{\link{dofuture}}, called in particular by its bootstrap procedures. Which one is used is determined by \code{spaMM.options(wrap_parallel="dopar")} (default) or \code{spaMM.options(wrap_parallel="dofuture")}. Depending on arguments, either serial computation (default), a socket cluster (parallelisation default), or a fork cluster (available in linux and alike operating systems) can be used. 

\code{dopar} is based on a patchwork of backends: for socket clusters, depending whether the \code{doSNOW} package is attached, \code{foreach} or \code{pbapply} is called (\code{doSNOW} allows more efficient load balancing than \code{pbapply}); for fork clusters, \code{parallel::mclapply} is used. By contrast, \code{dofuture} is based only on the \code{future} and \code{future.apply} packages, ensuring identical control of random number generator across these different cases, hence repeatable results across them. This does \bold{not} make a difference for bootstrap computations in spaMM as the bootstrap samples are never simulated in parallel: only refitting the models is performed in parallel, and fit results do not depend on random numbers. Further, the \code{future}-based code for socket clusters appears significantly slower than the one used by \code{dopar}. For these reasons, the latter function is used by default by spaMM.


}
