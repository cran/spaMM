\name{inner-vs-outer}
\alias{inner-vs-outer}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Inner iterations versus outer optimization
}
\description{
For estimation of \code{phi} and \code{lambda}, two approaches are implemented: generic optimization functions to maximize the log-likelihood (or restricted log-likelihood in REML case), and iterative algorithms described by Lee et al. (2006). This distinction is extended to random-coefficient parameters (\code{ranCoefs}) and to the correlation parameter of \code{\link{adjacency}} models in \pkg{spaMM}, but other correlation parameters of random effects are estimated only by generic optimization. Fixed-effect coefficients (\eqn{\beta}) are always estimated by iterative methods. The iterative method is well suited for estimating parameters of heteroscedastic (or \dQuote{structured-dispersion}) models, but it then requires computations which can be slow for large datasets. 

The \code{HLfit} and \code{HLCor} functions only use the iterative algorithms. Other fitting functions can call both approaches, internally using \code{HLfit} and \code{HLCor} to perform any iterative estimation. For this reason, all fitting functions, \bold{except} \code{HLfit} and \code{HLCor}, are called \bold{\dQuote{outer} optimization functions} in this documentation, possibly calling \bold{\dQuote{inner} iterative algorithms}. Likewise, parameters may be \bold{\dQuote{outer}}- or \bold{\dQuote{inner}}-estimated.

\code{fitme} by default tries to select the fastest method when both can be applied. Therefore, \code{fitme} tends to use generic generic optimization, rather than iterative algorithms, for large datasets, unless there is a non-trivial \code{resid.model}. \code{corrHLfit} tends to use iterative methods as much as possible.

Compared to generic optimization, the iterative algorithms provide more information, such as the \dQuote{cond. SE} of the dispersion parameter estimates (about which see warning in Details of \code{\link{HLfit}}). To force the evaluation of such information after an outer-optimization by an outer-optimizing function, use the \code{control$refit} argument (see Examples in \code{\link{fitme}} documentation). Alternatively (and possibly of limited use), one can force inner-optimization of \code{lambda} for a given random effect, or of \code{phi}, by setting it to \code{NaN} in \code{init} (see further Example using \sQuote{blackcap} data). The same syntax may be used for \code{phi}.

}
\details{
The results of REML fits of non-gaussian mixed models with both inner- and outer-estimated parameters may (generally slightly) differ according to which parameters are inner-estimated, when the fits include both some dispersion estimate(s) (say \eqn{\hat{\lambda}}),  which maximize restricted likelihood for given fixed effects, and the fixed effects estimates (\eqn{\hat{\beta}}, which maximize marginal likelihood for given \eqn{\hat{\lambda}}). Fit results then depend whether \eqn{\lambda} is inner- or outer-estimated, if \eqn{\hat{\lambda}} depends on \eqn{\hat{\beta}} and conversely. The latter condition may have been ignored on the grounds of asymptotic arguments, but is nevertheless observed. 

Explicit values of the \code{init} argument of the fitting functions may be used to control the use of inner iterations, and may be used to ensure better consistency of results over time if selection of methods changes changes over successive versions of \code{spaMM}: use \code{NaN} as noted above to enforce inner-optimization, and \code{NA} or a numeric value to enforce outer-optimization.  
}
\references{
  Lee, Y., Nelder, J. A. and Pawitan, Y. (2006). Generalized linear models with random effects: unified analysis via h-likelihood. Chapman & Hall: London.
}
\examples{
###  Contrasting different optimization methods:

## Toy data:
set.seed(123)
gr <- data.frame(y=rgamma(100,shape=9/2,scale=2/3))
# => Gamma deviates with mean mu=3 and variance=2, 
#  ie. phi= var/mu^2= 2/9 

## Fits

# Here fitme uses HLfit methods which provide cond. SE for phi by default:
fitme(y~1,data=gr,family=Gamma(log))

# To force outer optimization of phi, use the init argument:
fitme(y~1,data=gr,family=Gamma(log),init=list(phi=1))

# To obtain cond. SE for phi after outer optimization, use the 'refit' control:
fitme(y~1,data=gr,family=Gamma(log),,init=list(phi=1),
      control=list(refit=list(phi=TRUE))) ## or ...refit=TRUE...

### Outer-optimization is not necessarily the best way to find a global maximum, 
#  particularly when there is little statistical information in the data:  

if (spaMM.getOption("example_maxtime")>1.6) {
  data("blackcap")
  fitme(migStatus ~ means+ Matern(1|longitude+latitude),data=blackcap) # poor
  #  Compare with the following two ways of avoiding outer-optimization of lambda:
  corrHLfit(migStatus ~ means+ Matern(1|longitude+latitude),data=blackcap,
            method="ML")
  fitme(migStatus ~ means+ Matern(1|longitude+latitude),data=blackcap, 
        init=list(lambda=NaN))
}
}
