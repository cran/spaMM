\name{get_RLRsim_args}
\alias{get_RLRTSim_args}
\alias{get_RLRsim_args}
\title{
Extractors of arguments for functions from package RLRsim
}
\description{
\code{get_RLRsim_args} extracts a list of arguments suitable for a call to \code{RLRsim::RLRTSim()} or \code{RLRsim::LRTSim()}. These functions use an efficient simulation procedure to compute restricted or marginal likelihood ratio tests, respectively, comparing a fixed-effect model and a mixed-effect model with one random effect. They are notably used to test for the presence of one random effect, although the models compared by marginal likelihood (\code{LRTSim()}) may differ both in their random and in their fixed effects (as shown in the Example). The tests are exact for small samples (up to simulation error) for LMMs with no free parameters in the random effect (beyond the variance being tested), so not for time-series or spatial models with estimated correlation parameters. Heteroscedasticity of the residuals or of the random effect variance are also not taken into account by the simulation procedure (see Value field below for an hint why this is so). 

\code{get_RLRTSim_args} is the older extractor, originally for \code{RLRsim::RLRTSim()} only, now handling also ML fits with a warning. 
}
\usage{
get_RLRsim_args(object, verbose=TRUE, REML=NA, ...)
get_RLRTSim_args(object, verbose=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of class \code{HLfit}, as returned by the fitting functions in \code{spaMM}, for the more complete model to be compared.}
  \item{verbose}{NA or boolean; Whether to display some message or not.}
  \item{REML}{For programming purposes, not documented.}
  \item{\dots}{Additional arguments (currently not used).}
}
\value{
A list of arguments for a call to \code{RLRsim::RLRTSim()}  or \code{RLRsim::LRTSim()}. The main arguments are the design matrices of the more complete model (the one for the fixed effect, and the \bold{ZAL} matrix detailed in \code{\link{random-effects}}). The models handled by the testing procedure are the ones that are sufficiently characterized by these two matrices.
}
\note{
The inconsistent capitalisation of 's' in the function names is consistent with the inconsistencies in the \code{RLRsim} package.
}
\references{
Crainiceanu, C. and Ruppert, D. (2004) Likelihood ratio tests in
  linear mixed models with one variance component, \emph{Journal of the Royal
  Statistical Society: Series B},\bold{66},165--185.
}
\examples{
\dontrun{
## Derived from example in RLRsim::LRTSim
 set.seed(123)
 dat <- data.frame(g = rep(1:10, e = 10), x = (x<-rnorm(100)), 
                   y = 0.1 * x + rnorm(100))
 m <- fitme(y ~ x + (1|g), data=dat)
 m0 <- fitme(y ~ 1, data=dat) 
 (obs.LRT <- 2*(logLik(m)-logLik(m0)))
 args <- get_RLRsim_args(m)
 sim.LRT <- do.call(RLRSim::RLRTsim, args )
 (pval <- (sum(sim.LRT >= obs.LRT) + 1) / (length(sim.LRT) + 1))
}
}
\keyword{htest}
\keyword{regression}