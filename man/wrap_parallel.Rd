\name{wrap_parallel}
\alias{wrap_parallel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Selecting interfaces for parallelisation
}
\description{
spaMM implements three interfaces for parallelisation. Depending on their arguments, either serial computation (default), a socket cluster (parallelisation default), or a fork cluster (available in linux and alike operating systems) can be used by all interfaces. 

\code{\link{dopar}} is called by default by its bootstrap procedures, and \code{\link{dofuture}} has been developed as an alternative, whose use is controlled by \code{spaMM.options(wrap_parallel="dofuture")} (versus the default, \code{spaMM.options(wrap_parallel="dopar")}. \code{\link{combinepar}} is the third and more recent interface; it is not a formally supported \code{wrap_parallel} option because its additional functionalities are of no use in \pkg{spaMM}'s bootstrap procedures.

\code{dopar} is based on a patchwork of backends: for socket clusters, depending whether the \pkg{doSNOW} package is attached, \code{foreach} or \code{pbapply} is called (\pkg{doSNOW} allows more efficient load balancing than \code{pbapply}); for fork clusters, \code{parallel::mclapply} is used. This makes it impossible to ensure consistency of options accross computation environments, notably of enforcing the \code{.combine} control of \code{foreach}; and this makes it error-prone to ensure identical control of random number generators in all cases (although \code{dopar} and \code{combinepar} still aim to ensure the latter control).  

By contrast, \code{dofuture} is based only on the \pkg{future} and \pkg{future.apply} packages, in principle allowing a single syntax to control of random number generator across the different cases, hence repeatable results across them. This does \bold{not} make a difference for bootstrap computations in \pkg{spaMM} as the bootstrap samples are never simulated in parallel: only refitting the models is performed in parallel, and fit results do not depend on random numbers. Further, the \pkg{future}-based code for socket clusters appears significantly slower than the one used by \code{dopar}. For these reasons, the latter function is used by default by \pkg{spaMM}.

\code{combinepar} is a third and more recent approach designed to address the other issue: it always uses \pkg{foreach} so that the \code{.combine} control is consistently enforced. It uses \pkg{future} only when no alternative is available to produce a progress bar (namely, for socket clusters when \pkg{doSNOW} is not available).

}
