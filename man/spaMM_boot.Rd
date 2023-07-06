\name{spaMM_boot}
\alias{spaMM_boot}
\alias{spaMM2boot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parametric bootstrap
}
\description{
\code{spaMM_boot} simulates samples from a fit object inheriting from class \code{"HLfit"}, as produced by spaMM's fitting functions, and applies a given function to each simulated sample. Parallelization is supported (see Details). 

\code{spaMM2boot} is similar except that it assumes that the original model is refitted on the simulated data, and the given function is applied to the refitted model, and the value is in a format directly usable as input for \code{boot::boot.ci}. 

Both of these functions can be used to apply standard parametric bootstrap procedures. \code{spaMM_boot} is suitable for more diverse applications, e.g. to fit by one model some samples simulated under another model (see Example).
}
\usage{
spaMM_boot(object, simuland, nsim, nb_cores=NULL, seed=NULL,
           resp_testfn=NULL, control.foreach=list(),
           debug. = FALSE, type, fit_env=NULL, cluster_args=NULL,
           showpbar= eval(spaMM.getOption("barstyle")),
           boot_samples=NULL,
           ...)
spaMM2boot(object, statFUN, nsim, nb_cores=NULL, seed=NULL,
           resp_testfn=NULL, control.foreach=list(),
           debug. = FALSE, type="marginal", fit_env=NULL, 
           cluster_args=NULL, showpbar= eval(spaMM.getOption("barstyle")),
           boot_samples=NULL,
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
  \item{statFUN}{
The function to apply to each fit object for each simulated sample. See Details for requirements of this function.
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
  Boolean (or integer, interpreted as boolean). For debugging purposes, given that \code{spaMM_boot} does not stop when the fit of a bootstrap replicate fails. Subject to changes with no or little notice. In serial computation, \code{debug.=2} will stop on an error. In parallel computation, this would be ignored. The effect of \code{debug.=TRUE} depends on what \code{simuland} does of it. The default \code{simuland} for likelihood ratio testing functions, \code{\link{eval_replicate}}, shows how \code{debug.} can be used to control a call to \code{dump.frames} (however, debugging user-defined functions by such a call does not require control by \code{debug.}).  
}
\item{type}{
  Character: passed to \code{simulate.HLfit}. Defaults, with a warning, to \code{type="marginal"} in order to replicate the behaviour of previous versions of \code{spaMM_boot}. This is an appropriate default for various parametric bootstrpa analyses, but not necessarily the appropriate \code{type} for all possible uses. See Details of \code{\link{simulate.HLfit}} for other implemented options.
}  
\item{fit_env}{
  An environment or list containing variables necessary to evaluate \code{simuland} on each sample, and not included in the fit \code{object}. E.g., use \code{fit_env=list(phi_fix=phi_fix)} if the fit assumed \code{fixed=list(phi=phi_fix)}: the name in \code{list(phi_fix=<.>)} must be the name of the object that will be sought by the called process when interpreting \code{fixed=list(phi=phi_fix)} (if still unsure about the proper syntax, see the \code{\link[parallel]{clusterExport}} documentation, as \code{fit_env} is used in the following context: \code{parallel::clusterExport(cl=<cluster>, varlist=ls(fit_env), envir=fit_env))}.   
} 
\item{cluster_args}{
  \code{NULL} or a \code{list} of arguments, passed to \code{\link[parallel]{makeCluster}}. 
}  
  \item{showpbar}{
Controls display of progress bar. See \code{\link{barstyle}} option for details.
}
  \item{boot_samples}{
NULL, or a list as returned by \code{\link{simulate4boot}}. This may be used to provide precomputed bootstrap samples from the fitted model.
}
  \item{\dots}{
  Further arguments passed to the \code{simuland} function.
}
}
\details{
The \code{simuland} function must take as first argument a vector of response values, and may have other arguments including \sQuote{\ldots}. When required, these additional arguments must be passed through the \sQuote{\ldots} arguments of \code{spaMM_boot}. Variables needed to evaluate them must be available from within the \code{simuland} function or otherwise provided as elements of \code{fit_env}.

The \code{statFUN} function must take as first argument a fit object, and may have other arguments including \sQuote{\ldots} handled as for \code{simuland}.

\code{spaMM_boot} handles parallel backends with different features. \code{pbapply::pbapply} has a very simple interface (essentially equivalent to \code{apply}) and provides progress bars, but (in version 1.4.0, at least) does not have efficient load-balancing. \code{doSNOW} also provides a progress bar and allows more efficient load-balancing, but its requires \code{foreach}. \code{foreach} handles errors differently from \code{pbapply} (which will simply stop if fitting a model to a bootstrap replicate fails): see the \code{foreach} documentation.

\code{spaMM_boot} calls \code{simulate.HLfit} on the fit \code{object} and applies \code{simuland} on each column of the matrix returned by this call. 
\code{simulate.HLfit} uses the \code{type} argument, which must be explicitly provided.
}
\value{\code{spaMM_boot} returns a list, with two elements (unless \code{debug.} is \code{TRUE}): 
  \describe{
    \item{bootreps}{\code{nsim} return values in the format returned either by \code{apply} or \code{parallel::parApply} or by \code{foreach::`\%dopar\%`} as controlled by \code{control.foreach$.combine} (which is here \code{"rbind"} by default).}
    \item{RNGstate}{the state of \code{.Random.seed} at the beginning of the sample simulation}.  
  }
  
  \code{spaMM2boot} returns a list suitable for use by \code{boot.ci}, with elements:  
  \describe{
    \item{t}{\code{nsim} return values of the simulated statistic (in matrix format).}
    \item{t0}{\code{nsim} return the value of \code{statFUN} from the original fit.}
    \item{sim}{The simulation type (\code{"parametric"}).}
    \item{R}{\code{nsim}}
    \item{.Random.seed}{the state of \code{.Random.seed} at the beginning of the sample simulation}.  
  }
  (other elements of an object of class \code{\link{boot}} are currently not included.)

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
  spaMM_boot(lrt$nullfit, simuland = myfun, nsim=1, lrt=lrt, 
             type ="marginal")[["bootreps"]] 
  
  # Return only a model coefficient for each fit: 
  spaMM_boot(lrt$nullfit, simuland = myfun, nsim=7,
             what=quote(function(res) fixef(res)[2L]), 
             lrt=lrt, type ="marginal")[["bootreps"]]       
  
  \dontrun{
    # Parametric bootstrap by spaMM2boot() and spaMM_boot():
    boot.ci_info <- spaMM2boot(lrt$nullfit, statFUN = function(refit) fixef(refit)[1], 
                               nsim=99, type ="marginal")
    boot::boot.ci(boot.ci_info, , type=c("basic","perc","norm"))
    
    nullfit <- lrt$nullfit
    boot_t <- spaMM_boot(lrt$nullfit, simuland = function(y, nullfit) {
      refit <- update_resp(nullfit, y)
      fixef(refit)[1]
    }, nsim=99, type ="marginal", nullfit=nullfit)$bootreps
    boot::boot.ci(list(R = length(boot_t), sim="parametric"), t0=fixef(nullfit)[1], 
                  t= t(boot_t), type=c("basic","perc","norm"))


  }           
}
}
\keyword{ htest }
