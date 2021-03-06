\name{spaMM_boot}
\alias{spaMM_boot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parametric bootstrap
}
\description{
This simulates samples from a fit object inheriting from class \code{"HLfit"}, as produced by spaMM's fitting functions, and applies a given function to each simulated sample. Parallelization is supported (see Details). A typical usage of the parametric bootstrap is to fit by one model some samples simulated under another model (see Example).
}
\usage{
spaMM_boot(object, simuland, nsim, nb_cores=NULL, seed=NULL,
           resp_testfn=NULL, control.foreach=list(),
           debug. = FALSE, type, fit_env=NULL, cluster_args=NULL,
           showpbar= eval(spaMM.getOption("barstyle")),
           ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  The fit object to simulate from.
}
  \item{simuland}{
The function to apply to each simulated sample. See Details for requirements of this function.
}
  \item{nsim}{
Number of samples to simulate and analyze.
}
  \item{nb_cores}{
  Number of cores to use for parallel computation. The default is \code{spaMM.getOption("nb_cores")}, and 1 if the latter is NULL. \code{nb_cores=1} prevents the use of parallelisation procedures.
}
\item{seed}{Passed to \code{\link{simulate.HLfit}}}
  \item{resp_testfn}{Passed to \code{simulate.HLfit}; NULL, or a function that tests a condition which simulated samples should satisfy. This function takes a response vector as argument and return a boolean (TRUE indicating that the sample satisfies the condition).     
  }
  \item{control.foreach}{list of control arguments for \code{foreach}. These include in particular \code{.combine} (with default value \code{"rbind"}), and \code{.errorhandling} (with default value \code{"remove"}, but \code{"pass"} is quite useful for debugging).}
  \item{debug.}{
  Boolean (or integer, interpreted as boolean). For debugging purposes, particularly from parallel computations. The effect of \code{debug.=TRUE} depends on what \code{simuland} does of it. The default \code{simuland} for likelihood ratio testing functions, \code{\link{eval_replicate}}, shows how \code{debug.} can be used.
}
\item{type}{
  Character: passed to \code{simulate.HLfit}. Defaults, with a warning, to code{type="marginal"} in order to replicate the behaviour of previous versions. But this is not necessarily the appropriate \code{type} for all possible uses. See Details of \code{\link{simulate.HLfit}} for other implemented options.
}  
\item{fit_env}{
  An environment or list containing variables necessary to evaluate \code{simuland} on each sample, and not included in the fit \code{object}. E.g., use \code{fit_env=list(phi_fix=phi_fix)} if the fit assumed \code{fixed=list(phi=phi_fix)} 
} 
\item{cluster_args}{
  \code{NULL} or a \code{list} of arguments, passed to \code{\link[parallel]{makeCluster}}. 
}  
  \item{showpbar}{
Controls display of progress bar. See \code{\link{barstyle}} option for details.
}
  \item{\dots}{
  Further arguments passed to the \code{simuland} function.
}
}
\details{
\code{spaMM_boot} handles parallel backends with different features. \code{pbapply::pbapply} has a very simple interface (essentially equivalent to \code{apply}) and provides progress bars, but (in version 1.4.0, at least) does not have efficient load-balancing. \code{doSNOW} also provides a progress bar and allows more efficient load-balancing, but its requires \code{foreach}. \code{foreach} handles errors diferently from \code{pbapply} (which will simply stop if fitting a model to a bootstrap replicate fails): see the \code{foreach} documentation.

\code{spaMM_boot} calls \code{simulate.HLfit} on the fit \code{object} and applies \code{simuland} on each column of the matrix returned by this call. 
\code{simulate.HLfit} uses the \code{type} argument, which must be explicitly provided.

The \code{simuland} function must take as first argument a vector of response values, and may have other arguments including \sQuote{\ldots}. When required, these additional arguments must be passed through the \sQuote{\ldots} arguments of \code{spaMM_boot}. Variables needed to evaluate them must be available from within the \code{simuland} function or otherwise provided as elements of \code{fit_env}.
}
\value{A list, with two elements (unless \code{debug.} is \code{TRUE}): 
  \describe{
    \item{bootreps}{\code{nsim} return values in the format returned either by \code{apply} or \code{parallel::parApply} or by \code{foreach::`\%dopar\%`} as controlled by \code{control.foreach$.combine} (which is here \code{"rbind"} by default).}
    \item{RNGstate}{the state of \code{.Random.seed} at the beginning of the sample simulation}.  
  }
}
\examples{
if (spaMM.getOption("example_maxtime")>7) {
 data("blackcap")
 
 # Generate fits of null and full models:
 lrt <- fixedLRT(null.formula=migStatus ~ 1 + Matern(1|longitude+latitude),
                 formula=migStatus ~ means + Matern(1|longitude+latitude), 
                 method='ML',data=blackcap)

 # The 'simuland' argument: 
 myfun <- function(y, what=NULL, lrt, ...) { 
    data <- lrt$fullfit$data
    data$migStatus <- y ## replaces original response (! more complicated for binomial fits)
    full_call <- getCall(lrt$fullfit) ## call for full fit
    full_call$data <- data
    res <- eval(full_call) ## fits the full model on the simulated response
    if (!is.null(what)) res <- eval(what)(res=res) ## post-process the fit
    return(res) ## the fit, or anything produced by evaluating 'what'
  }
  # where the 'what' argument (not required) of myfun() allows one to control 
  # what the function returns without redefining the function.
  
  # Call myfun() with no 'what' argument: returns a list of fits 
  spaMM_boot(lrt$nullfit, simuland = myfun, nsim=1, lrt=lrt, type ="marginal")[["bootreps"]] 
  
  # Return only a model coefficient for each fit: 
  spaMM_boot(lrt$nullfit, simuland = myfun, nsim=7,
               what=quote(function(res) fixef(res)[2L]), lrt=lrt, type ="marginal")[["bootreps"]]       
}
}
\keyword{ htest }
