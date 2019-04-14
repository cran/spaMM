\name{spaMM_boot}
\alias{spaMM_boot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parametric bootstrap
}
\description{
This simulates samples from a fit object inheriting from class \code{"HLfit"}, as produced by spaMM's fitting function, and applies a given function to each simulated sample. Parallelization is supported (see Details). A typical usage of the parametric bootstrap is to fit by one model samples produced under another model (see Example).
}
\usage{
spaMM_boot(object, simuland, nsim, nb_cores = NULL, resp_testfn=NULL, 
           control.foreach=list(), debug. = FALSE, type, ...)
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
  \item{resp_testfn}{Passed to \code{simulate.HLfit}; NULL, or a function that tests a condition which simulated samples should satisfy. This function takes a response vector as argument and return a boolean (TRUE indicating that the sample satisfies the condition).     
  }
  \item{control.foreach}{list of control arguments for \code{foreach}. These include in particular \code{.combine} (with default value \code{"rbind"}), and \code{.errorhandling} (with default value \code{"remove"}, but \code{"pass"} is quite useful for debugging).}
  \item{debug.}{
  Boolean; only for debugging purposes, particularly in parallel runs, where \code{debug.=TRUE} allows useful debugging info to be returned. If \code{debug.=FALSE}, the returned \code{bootreps} will contain \code{NA} for fits that fail, and if \code{debug.=2} an error will be thrown for fits that fail (only useful in serial computations). 
}
\item{type}{
  For development purposes, not further documented. % F I X M E maybe remove
}  
  \item{\dots}{
  Further arguments passed (in a non-standard way) to the \code{simuland} function.
}
}
\details{
\code{spaMM_boot} handles parallel backends with different features. \code{pbapply::pbapply} has a very simple interface (essentially equivalent to \code{apply}) and provides progress bars, but (currently: version 1.4.0) does not have efficient load-balancing. \code{doSNOW} also provides a progress bar and allows more efficient load-balancing, but its requires \code{foreach}, whose handling of '\dots' arguments is tortuous. \code{foreach} will be used if \code{doSNOW} is loaded; then, some of the '\dots' arguments may need to be quoted (see Example). \code{foreach} also handles errors diferently from \code{pbapply} (which will simply stop if fitting a model to a bootstrap replicate fails): see the \code{foreach} documentation.

\code{spaMM_boot} calls \code{simulate.HLfit} on the fit \code{object} and applies \code{simuland} on each column of the matrix returned by this call. The \code{simuland} function must take as first argument a vector of response values, and may use \ldots to pass additional arguments. 
The default \code{simuland} function is \code{.eval_replicate()}, and an alternative function \code{.eval_replicate2()} is also provided. The latter function is slower, as it refits the models compared with different initial values for random-effect parameters, which is useful in some difficult cases where initial values matter. 

Advanced users can define their own \code{simuland} function. An example is provided in the file \code{tests/testthat/test-LRT-boot.R}, using \ldots to pass additional arguments beyond response values. Alternatively, \code{.eval_replicate()} can be used as a template. It has no \ldots argument, as essential arguments are passed through its environment. \code{spaMM_boot} redefines the environment of any \code{simuland} for that purpose. This implies that users but should not assume that they can control their own \code{simuland} function's environment (except its \code{isdebugged()} status). 
}
\value{A list, with two elements (unless \code{debug.} is \code{TRUE}): 
  \itemize{
    \item \code{bootreps}, \code{nsim} return values in the format returned either by \code{apply} or \code{parallel::parApply} or by \code{foreach::`\%dopar\%`} as controlled by \code{control.foreach$.combine}. If \code{simuland} returns a vector, \code{spaMM_boot} should effectively \code{rbind} the results by default, returning an \code{nsim}-row matrix in all cases. From spaMM 2.5.6, if \code{simuland} returns a 1-row data frame, \code{spaMM_boot} \code{rbind}s the results into a \code{nsim}-row data frame in all cases. The results may not be consistent among parallel backends in other cases, and may change in later versions, so users should stick to one of these two cases as much as possible.
    \item \code{RNGstate}, the state of \code{.Random.seed} at the beginning of the simulation.  
  }
}
\examples{
if (spaMM.getOption("example_maxtime")>10) {
 data("blackcap")
 
 # Generate fits of null and full models:
 lrt <- fixedLRT(null.formula=migStatus ~ 1 + Matern(1|latitude+longitude),
       formula=migStatus ~ means + Matern(1|latitude+longitude), 
       HLmethod='ML',data=blackcap)

 # The 'simuland' argument: 
 myfun <- function(y, what=NULL, lrt, ...) { 
    data <- lrt$fullfit$data
    data$migStatus <- y ## replaces original response (! more complicated for binomial fits)
    full_call <- getCall(lrt$fullfit) ## call for full fit
    full_call$data <- data
    res <- eval(full_call) ## fits the full model on the simulated response
    if (!is.null(what)) res <- eval(what) ## post-process the fit
    return(res) ## the fit, or anything produced by evaluating 'what'
  }
  # where the 'what' argument (not required) of myfun() allows one to control 
  # what the function returns without redefining the function.
  
  # Call myfun() with no 'what' argument: returns a list of fits 
  spaMM_boot(lrt$nullfit, simuland = myfun, nsim=1, lrt=lrt)[["bootreps"]] 
  
  # Return only a model coefficient for each fit: 
  spaMM_boot(lrt$nullfit, simuland = myfun, nsim=7,
               what=quote(fixef(res)[2L]), lrt=lrt)[["bootreps"]]       
}
}
\keyword{ htest }
