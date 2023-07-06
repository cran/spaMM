\name{simulate.HLfit}
\alias{simulate.HLfit}
\alias{simulate.HLfitlist}
\alias{simulate_ranef}
\alias{simulate}
\alias{simulate4boot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Simulate realizations of a fitted model.
}
\description{
  From an HLfit object, \code{simulate.HLfit} function generates new samples given the estimated fixed effects 
  and dispersion parameters. Simulation may be unconditional (the default, useful in many applications of parametric bootstrap), or conditional on the predicted values of random effects, or may draw from the conditional distribution of random effects given the observed response.
  Simulations may be run for the original values of fixed-effect predictor variables and of random effect levels (spatial locations for spatial random effects), or for new values of these. 
  
  \code{simulate4boot} is a wrapper around \code{simulate.HLfit} that can be used to precompute the bootstrap samples to be used by \code{\link{spaMM_boot}} or \code{\link{spaMM2boot}} through their \code{boot_samples} argument (and is called internally by these functions when \code{boot_samples} is NULL).
  
  \code{simulate_ranef} will only simulate and return a vector of random effects, more specifically some elements of \bold{b} in the standard form \code{offset}+ \bold{X}\eqn{\beta} + \bold{Z b} for the linear predictor.
}
\usage{
\method{simulate}{HLfit}(object, nsim = 1, seed = NULL, newdata = NULL, 
                         type = "marginal", re.form, conditional = NULL, 
                         verbose = c(type=TRUE,
                                     showpbar=eval(spaMM.getOption("barstyle"))), 
                         sizes = NULL, resp_testfn = NULL, phi_type = "predict", 
                         prior.weights = object$prior.weights, variances=list(), ...)

\method{simulate}{HLfitlist}(object, nsim = 1, seed = NULL, 
                             newdata = object[[1]]$data, sizes = NULL, ...)

simulate4boot(object, nsim, seed=NULL, resp_testfn=NULL, type="marginal", 
              showpbar=eval(spaMM.getOption("barstyle")))
simulate_ranef(object, which=NULL, newdata = NULL, nsim = 1L)
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
  A data frame closely matching the original data, except that response values are not needed. May provide new values of fixed predictor variables, new spatial locations, new individuals within a block, or new values of the LHS in random-effect terms of the form \code{(<LHS>|<RHS>)}.   
}
\item{re.form}{ formula for random effects to condition on. Default behaviour depends on the \code{type} argument. The joint default is the latter's default, i.e., unconditional simulation. \code{re.form} is currently ignored when 
\code{type="predVar"} %also (ranef|response) but this type is not mentioned in the doc
(with a warning). Otherwise, \code{re.form=NULL} conditions on all random effects (as \code{type="residual"} does), and \code{re.form=NA} conditions on none of the random effects (as \code{type="marginal"} or \code{re.form=~0} do).  
}
  \item{type}{
    character string specifying which uncertainties are taken into account in the linear predictor and notably in the random effect terms. Whatever the \code{type}, the residual variance is always accounted in the simulation output. \code{"marginal"} accounts for the marginal variance of the random effect (and, by default, also for the uncertainty in fixed effects); \code{"predVar"} accounts for the conditional distribution of the random effects given the data (see Details); and \code{"residual"} should perhaps be \code{"none"} as no uncertainty is accounted in the linear predictor: the simulation variance is only the residual variance of the fitted model.
  }
  \item{conditional}{
    Obsolete and will be deprecated. Boolean; TRUE and FALSE are equivalent to \code{type="residual"} and \code{type="marginal"}, respectively. 
  }
  \item{verbose}{
    Either a single boolean (which determines \code{verbose[["type"]]}, or a vector of booleans with possible elements \code{"type"} (to display basic information about the type of simulation) and \code{"showpbar"} (see \code{predict(.,verbose)}). 
  }
  \item{sizes}{
  A vector of sample sizes to simulate in the case of a binomial fit. Defaults to the sizes in the original data. 
}
  \item{resp_testfn}{
  NULL, or a function that tests a condition which simulated samples should satisfy. This function takes a response vector as argument and return a boolean (TRUE indicating that the sample satisfies the condition).     
  }
  \item{phi_type}{
    Character string, either \code{"predict"} or one of the values possible for \code{type}. This controls the residual variance parameter \eqn{\phi}. The default is to use predicted \eqn{\phi} values from the fit, which are the fitted \eqn{\phi} values except when a structured-dispersion model is involved together with non-NULL \code{newdata}. However, when a structured-dispersion model is involved, it is also possible to simulate new \eqn{\phi} values, and for a mixed-effects structured-dispersion model, the same types of simulation controlled by \code{type} for the mean response can be performed as controlled by \code{phi_type}. For a fixed-effects structured-dispersion model, these types cannot be distinguished, and any \code{phi_type} distinct from \code{"predict"} will imply simulation under the fixed-effect model (see Examples). 
  }
  \item{prior.weights}{
  Prior weights that may be substituted to those of the original fit, with the same effect on the residual variance. 
}
\item{variances}{
  Used when \code{type="predVar"}: see Details. 
}  
\item{\dots}{
  further arguments passed to or from other methods; currently only passed to \code{predict} in a speculative bit of code (see Details). 
}
\item{which}{
  Integer, or integer vector: the random effect(s) (indexed as ordered as in the model formula) to be simulated. If NULL, all of them are simulated. 
}  
\item{showpbar}{
Controls display of progress bar. See \code{\link{barstyle}} option for details.
}  
}
\details{
\code{type="predVar"} accounts for the uncertainty of the linear predictor, by drawing new values of the predictor in a multivariate gaussian distribution with mean and covariance matrix of prediction. In this case, the user has to provide a \code{variances} argument, passed to \code{predict}, which controls what goes into this covariance matrix. For example, with \code{variances=list(linPred=TRUE,disp=TRUE)}), the covariance matrix takes into account the joint uncertainty in the fixed-effect coefficients and of any random effects given the response and the point estimates of dispersion and correlation parameters (\code{"linPred"} variance component), and in addition accounts for uncertainty in the dispersion parameters (effect of \code{"disp"} variance component as further described in \code{\link{predict.HLfit}}).  The total simulation variance is then the response variance. Uncertainty in correlation parameters (such a parameters of the Matern family) is not taken into account. The \code{"linPred"} uncertainty is known exactly in LMMs, and otherwise approximated as a Gaussian distribution with mean vector and covariance matrix given as per the Laplace approximation. 

\code{type="(ranef|response)"} can be viewed as a special version of \code{type="predVar"} where\cr 
\code{variances=list(linPred=TRUE,disp=FALSE)}) and only the uncertainty in the random effects is taken into account.

A full discussion of the merits of the different \code{type}s is beyond the scope of this documentation, but these different types may not all be useful. \code{type="marginal"} is typically used for computation of confidence intervals by parametric bootstrap methods. \code{type="residual"} is used by \code{\link{get_cPredVar}} for its evaluation of a bias term. The other \code{type}s may be used to simulate the uncertainty in the random effects, conditionally on the data, and may therefore be more akin to the computation of prediction intervals conditionally on an (unknown but inferred) realization of the random effects. But these should presumably not be used in a bootstrap computation of such intervals, as this would represent a double accounting of the uncertainty that the boostrap aims to quantify.
}
\value{
  \code{simulate.HLfit} returns a vector (if nsim=1) or a matrix with \code{nsim} columns, each containing simulated responses (or simulated random effects, for \code{simulated_ranef()}).
  
  \code{simulate4boot} returns a list with elements 
    \describe{
    \item{bootreps}{the result of \code{simulate.HLfit} as a matrix with \code{nsim} columns;}
    \item{RNGstate}{the state of \code{.Random.seed} at the beginning of the sample simulation.}  
  }

  The \code{simulate.HLfitlist} method returns a list of simulated responses.
}
\examples{
data("Loaloa")
HLC <- HLCor(cbind(npos,ntot-npos)~Matern(1|longitude+latitude),
           data=Loaloa,family=binomial(),
           ranPars=list(lambda=1,nu=0.5,rho=1/0.7)) 
simulate(HLC,nsim=2)

## Structured dispersion model 
data("wafers")
hl <- HLfit(y ~X1+X2+X1*X3+X2*X3+I(X2^2)+(1|batch),family=Gamma(log),
            resid.model = ~ X3+I(X3^2) ,data=wafers)
simulate(hl,type="marginal",phi_type="simulate",nsim=2)
simulate_ranef(hl,nsim=2)
}
\keyword{ datagen }
