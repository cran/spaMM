\name{simulate.HLfit}
\alias{simulate.HLfit}
\alias{simulate.HLfitlist}
\alias{simulate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Simulate realizations of a fitted model.
}
\description{
  From an HLfit object, \code{simulate.HLfit} function generates new samples given the estimated fixed effects 
  and dispersion parameters. Simulation may be conditional on the predicted values of random effects, or unconditional.
  Simulations may be run for the original values of fixed-effect predictor variables and of random effect levels (spatial locations for spatial random effects), or for new values of these. However, conditional simulation in mixed models with new such values is not yet handled.  

  %This does not yet work for new locations in a mixed model involving both spatial and non-spatial random effects. 
}
\usage{
\method{simulate}{HLfit}(object, nsim = 1, seed = NULL, newdata=NULL,
                         conditional=FALSE, verbose=TRUE,
                         sizes = NULL, ...)
\method{simulate}{HLfitlist}(object, nsim = 1, seed = NULL, 
                             newdata=object[[1]]$data, sizes=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  The return object of HLfit or similar function.
}
  \item{nsim}{
 number of response vectors to simulate.  Defaults to '1'.
}
  \item{seed}{
  A seed for \code{\link{set.seed}}. If such a value is provided, the initial state of the random number generator at a global level is restored on exit from simulate.
}
  \item{newdata}{
  A data frame closely matching the original data, except that response values are not needed. May provide new values of fixed predictor variables, new spatial locations, or new individuals within a block.   
}
  \item{conditional}{
    Boolean; whether to run conditional or unconditional simulation. 
  }
  \item{verbose}{
    Boolean; whether to print some information or not. 
  }
  \item{sizes}{
  A vector of sample sizes to simulate in the case of a binomial fit. Defaults to the sizes in the original data. 
}
  \item{\dots}{
  further arguments passed to or from other methods. 
}
}
\value{
  For the \code{HLfitlist} method (i.e., the result of a multinomial fit), a list of simulated responses.
  Otherwise, a vector (if nsim=1) or a matrix with nsim columns, each containing a simulated response.
}
\examples{
data(Loaloa)
HLC <- HLCor(cbind(npos,ntot-npos)~Matern(1|longitude+latitude),
           data=Loaloa,family=binomial(),
           ranPars=list(lambda=1,nu=0.5,rho=1/0.7)) 
simulate(HLC,nsim=2)
}
\keyword{ datagen }
