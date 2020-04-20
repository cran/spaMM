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
  and dispersion parameters. Simulation may be unconditional (the default, useful in many applications of parametric bootstrap), or conditional on the predicted values of random effects, or may draw from the conditional distribution of random effects given the observed response.
  Simulations may be run for the original values of fixed-effect predictor variables and of random effect levels (spatial locations for spatial random effects), or for new values of these. 
}
\usage{
\method{simulate}{HLfit}(object, nsim = 1, seed = NULL, newdata = NULL, 
                         type = "marginal", re.form, conditional = NULL, 
                         verbose = c(type=TRUE,showpbar=interactive()), 
                         sizes = NULL, resp_testfn = NULL, phi_type = "predict", 
                         prior.weights = object$prior.weights, variances=list(), ...)
\method{simulate}{HLfitlist}(object, nsim = 1, seed = NULL, 
                             newdata = object[[1]]$data, sizes = NULL, ...)
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
\item{re.form}{ formula for random effects to condition on. Default behaviour depends on the \code{type} argument. The joint default is the latter's default, i.e., unconditional simulation. \code{re.form} is currently ignored when \code{type="Vlinpred"} (with a warning). Otherwise, \code{re.form=NULL} conditions on all random effects (as \code{type="residual"} does), and \code{re.form=NA} conditions on none of the random effects (as \code{type="marginal"} or \code{re.form=~0} do).  
}
  \item{type}{
    character string specifying the type of simulation for mixed models. \code{"marginal"} is for simulation from the marginal distribution of the random effect; \code{"residual"} accounts only for the residual variation of the fitted model; and the more speculative \code{"predVar"} accounts both for residual variation and for the uncertainty of the linear predictor (see Details). 
  }
  \item{conditional}{
    Obsolete and will be deprecated. Boolean; TRUE and FALSE are equivalent to \code{type="residual"} and \code{type="marginal"}, respectively. 
  }
  \item{verbose}{
    Either a single boolean (which determines \code{verbose[["type"]]}, or a vector of booleans with possible elements \code{"type"} (to display basic information about the type of simulation) and \code{"showpbar"} (see code{predict(.,verbose)}). 
  }
  \item{sizes}{
  A vector of sample sizes to simulate in the case of a binomial fit. Defaults to the sizes in the original data. 
}
  \item{resp_testfn}{
  NULL, or a function that tests a condition which simulated samples should satisfy. This function takes a response vector as argument and return a boolean (TRUE indicating that the sampel satisfies the condition).     
  }
  \item{phi_type}{
    Character string, either \code{"predict"} or one of the values possible for \code{type}. This controls the residual variance parameter \eqn{\phi}. The default is to use predicted \eqn{\phi} values from the fit, which are the fitted \eqn{\phi} values except when a structured-dispersion model is involved together with non-NULL \code{newdata}. However, when a structured-dispersion model is involved, it is also possible to simulate new \eqn{\phi} values, and for a mixed-effects structured-dispersion model, the same types of simulation controlled by \code{type} for the main response can be performed as controlled by \code{phi_type}. For a fixed-effects structured-dispersion model, these types cannot be distinguished, and any \code{phi_type} distinct from \code{"predict"} will imply simulation under the fixed-effect model (see Examples). 
  }
  \item{prior.weights}{
  Prior weights that may be substituted to those of the original fit, with the same effect on the residual variance. 
}
  \item{variances}{
  Used when \code{type="predVar"}: see Details. 
}  \item{\dots}{
  further arguments passed to or from other methods; currently only passed to \code{predict} in a speculative bit of code (see Details). 
}
}
\details{
\code{type="predVar"} accounts for the uncertainty of the linear predictor, by drawing new values of the predictor in a multivariate gaussian distribution with mean and covariance matrix of prediction. In this case, the user has to provide a \code{variances} argument, passed to \code{predict}. The covariance matrix represent the joint uncertainty in the fixed-effect coefficients and of any random effects given the response and the point estimates of dispersion and correlation parameters. This distribution is known exactly in LMMs, and otherwise approximated as a Gaussian distribution with mean vector and covariance matrix given as per the Laplace approximation. 

\code{type="(ranef|response)"} can be viewed as a special version of \code{type="predVar"} where \code{variances=list(linPred=TRUE,disp=FALSE)}) and only the uncertainty in the random effects is taken into account.

A full discussion of the merits of the different \code{type}s is beyond the scope of this documentation, but these different types may not all be useful. \code{type="marginal"} is typically used for computation of confidence intervals by parametric boostrap methods. \code{type="residual"} is used by \code{\link{get_cPredVar}} for its evaluation of a bias term. The other \code{type}s may be used to simulate the uncertainty in the random effects, conditionally on the data, and may thereofre be more akin to the computation of prediction intervals conditionnally on an (unknown but inferred) realization of the random effects. But these should presumably not be used in a bootstrap computation of such intervals, as this would represent a double accounting of the uncertainty that the boostrap aims to quantify.
}
\value{
  For the \code{HLfitlist} method (i.e., the result of a multinomial fit), a list of simulated responses.
  Otherwise, a vector (if nsim=1) or a matrix with nsim columns, each containing a simulated response.
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
}
\keyword{ datagen }
