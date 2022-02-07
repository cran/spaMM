\name{inits}
\alias{inits}
\title{
Controlling optimization strategy through initial values
}
\description{
Several parameters (notably the dispersion parameters: the variance of random effects and the residual variance parameter, if any) can be estimated either by iterative algorithms, or by generic optimization methods. The development of the \code{fitme} function aims to provide full control of the selection of algorithms. For example, if two random effects are fitted, then \code{init=list(lambda=c(NA,NaN))} enforces generic optimization for the first variance and iterative algorithms for the second. \code{init=list(lambda=c(0.1,NaN))} has the same effect and additionally provides control of the initial value for optimization (whereas \code{init.HLfit=list(lambda=c(NA,0.1))} will provide control of the initial value for iterations).

How to know which algorithm has been selected for each parameter? \code{fitme(., verbose=c(TRACE=TRUE))} shows successive values of the variables estimated by optimization (See Examples; if no value appears, then all are estimated by iterative methods). The first lines of the summary of a fit object should tell which variances are estimated by the \dQuote{outer} method.

\code{corrHLfit}, which uses inner optimization by default, can be forced to perform outer optimization. Its control is more limited, as \code{NA}s and \code{NaN}s are not allowed. Instead, only numeric values as in \code{init=list(lambda=0.1)} are allowed. 
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
# (but then it may be better to enforce it also for phi: init=list(lambda=NaN, phi=NaN))
#
# Use init to enforce generic optimization for lambda estimation,
#   and control initial value:      
fitme(passengers ~ month * year_z + AR1(1|time), data = air, 
      verbose=c(TRACE=TRUE), init=list(lambda=0.1))
      
# See help("multinomial") for more examples of control by initial values.       
}
}
