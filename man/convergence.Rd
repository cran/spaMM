\name{convergence}
\alias{convergence}
\title{
Assessing convergence for fitted models
}
\description{
spaMM fits may produce convergence warnings coming from \code{.check_conv_glm_reinit()}. These can generally be ignored (particularly when they show a small criterion, <1e-5).

spaMM fits may also produce messages pointing to slow convergence and drawing users here. These do not necessarily mean the fit is incorrect. Rather, they suggest that another fitting strategy could be tried. Keep in mind that several parameters (notably the dispersion parameters: the variance of random effects and the residual variance parameter, if any) can be estimated either by the iterative algorithms, or by generic optimization methods. In my experience, slow convergence happens in certain cases where a large random-effect variance is considered by the algorithm used. 
The development of the \code{fitme} function aims to provide full control of the selection of algorithms. For example, if two random effects are fitted, then \code{init=list(lambda=c(NA,NaN))} enforces generic optimization for the first variance and iterative algorithms for the second. \code{init=list(lambda=c(0.1,NaN))} has the same effect and additionnally provides control of the initial value for optimization (whereas \code{init.HLfit=list(lambda=c(NA,0.1))} will provide control of the initial value for iterations).

If the iterative algorithm is being used, then it is worth trying to use the generic optimization methods. In particular, if you used \code{HLfit}, try using \code{fitme}; if you already use \code{fitme}, try to enforce optimization of the random-effect variance(s). Conversely, if generic optimization is being used, the maximum lambda value could be controlled (say, \code{upper=list(lambda=c(10,NA))}), or the iterative algorithm can be called.

How to know which algorithm has been selected for each parameter? \code{fitme(., verbose=c(verbose=TRUE))} shows sucessive values of the variables estimated by optimization (See Examples; if no value appears, then all are estimated by iterative methods). 
}
\examples{
\dontrun{
air <- data.frame(passengers = as.numeric(AirPassengers),
                  year_z = scale(rep(1949:1960, each = 12)),
                  month = factor(rep(1:12, 12)))
air$time <- 1:nrow(air)
# Use verbose to find that lambda is estimated by optimization
fitme(passengers ~ month * year_z + AR1(1|time), data = air, 
      verbose=c(TRACE=TRUE)) 
# Use init to enforce iterative algorithm for lambda estimation:      
fitme(passengers ~ month * year_z + AR1(1|time), data = air, 
      verbose=c(TRACE=TRUE), init=list(lambda=NaN))
# (but then it may be better to enforce it also fro phi: init=list(lambda=NaN, phi=NaN))
#
# Use init to enforce generic optimization for lambda estimation,
#   and control initial value:      
fitme(passengers ~ month * year_z + AR1(1|time), data = air, 
      verbose=c(TRACE=TRUE), init=list(lambda=0.1))
}
}
