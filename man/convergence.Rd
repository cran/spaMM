\name{convergence}
\alias{convergence}
\title{
Assessing convergence for fitted models
}
\description{
spaMM fits can produce various convergence warnings or messages. 

Messages referring to convergence issues in initialization can generally be ignored but may help to diagnose other apparent problems, if any. 

Warnings referring to \code{.calc_dispGammaGLM} (for residual-dispersion fits) can generally be ignored when they show a small criterion (say <1e-5) but may otherwise suggest that the final fit did not optimize its objective. 

Messages pointing to slow convergence and drawing users to this doscumentation do not necessarily mean the fit is incorrect. Rather, they suggest that another fitting strategy could be tried. Keep in mind that several parameters (notably the dispersion parameters: the variance of random effects and the residual variance parameter, if any) can be estimated either by the iterative algorithms, or by generic optimization methods. In my experience, slow convergence happens in certain cases where a large random-effect variance is considered by the algorithm used. 

How to know which algorithm has been selected for each parameter? \code{fitme(., verbose=c(TRACE=TRUE))} shows successive values of the variables estimated by optimization (See Examples; if no value appears, then all are estimated by iterative methods). The first lines of the summary of a fit object should tell which variances are estimated by the \dQuote{outer} method.

If the iterative algorithm is being used, then it is worth trying to use the generic optimization methods. In particular, if you used \code{HLfit}, try using \code{fitme}; if you already use \code{fitme}, try to enforce generic optimization of the random-effect variance(s) (see \code{\link{inits}}). Conversely, if generic optimization is being used, the maximum lambda value could be controlled (say, \code{upper=list(lambda=c(10,NA))}), or the iterative algorithm can be called  (see \code{\link{inits}} again).

For the largest datasets, it may be worth comparing the speed of the \code{"spcorr"} and  \code{"spprec"} choices of the \code{\link{algebra}} control, in case spaMM has not selected the most appropriate by default. However, this will not be useful for geostatistical models with many spatial locations.
}
\examples{
# See help("inits") for examples of control by initial values.       
}
