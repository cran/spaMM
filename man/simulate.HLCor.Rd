\name{simulate.HLfit}
\alias{simulate.HLfit}
\alias{simulate.HLCor}
\alias{simulate.corrHLfit}
\alias{simulate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Simulate realizations of a fitted correlated mixed model.
}
\description{
  From an HLfit object, \code{simulate.HLfit} function generates new samples given the estimated fixed effects 
  and dispersion parameters, and given the design matrix for random effects (hence, given the correlation parameters).
  \code{simulate.HLCor} and \code{simulate.corrHLfit} are wrapper functions that extract from their respective input objects the HLfit object 
  and apply \code{simulate.HLfit} to it.  
}
\usage{
\method{simulate}{HLfit}(object, nsim = 1, seed = NULL, ...)
\method{simulate}{HLCor}(object,...)
\method{simulate}{corrHLfit}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  The return object of an HLfit, HLCor or corrHLfit call.
}
  \item{nsim}{
 number of response vectors to simulate.  Defaults to '1'.
}
  \item{seed}{
  A seed for \code{\link{set.seed}}. If such a value is provided, the initial state of the random number generator at a global level
  is restored on exit from simulate.
}
  \item{\dots}{
  further arguments passed to or from other methods. 
}
}
\value{
  A vector (if nsim=1) or a matrix with nsim columns, each containing a simulated response.
}
\examples{
data(Loaloa)
HLC <- HLCor(cbind(npos,ntot-npos)~Matern(1|longitude+latitude),
           data=Loaloa,family=binomial(),
           ranPars=list(lambda=1,nu=0.5,rho=1/0.7)) 
simulate(HLC,nsim=2)
}
\keyword{ distribution }
