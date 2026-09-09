\name{tweedie}
\alias{tweedie}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Tweedie (compound Poisson Gamma) family 
}
\description{
Returns a GLM \code{\link[stats]{family}} object for the Tweedie family of distributions (e.g., Dunn & Smyth, 2018). The variance of such distributions is related to the mean \eqn{\mu}: variance=\eqn{\phi\mu^p}, where \eqn{\phi} is a dispersion parameter (as in Gaussian and Gamma families), and \eqn{p} is the \code{index} parameter of the distribution, which can be specified as the \code{index} argument of the function, and is referred to in other contexts (such as initial values for optimization) as \code{Tw_index}. The index need or need not be specified when calling \code{tweedie}, depending on usage. The current implementation handles values \eqn{1<p<2}, for which the distributions are a Poisson sum of Gamma distributions.

The family link is controlled by another argument, \code{link}, that may specify a log link, or the exponent \eqn{q} of the link relationship \eqn{\eta=\mu^q}.

\code{\link[stats]{glm}} can be used with this family object and given \code{index}.

Fitting Tweedie models by \pkg{spaMM}'s fitting functions requires some functions from the \pkg{tweedieDistr} package, or else from the \pkg{tweedie} package, so either one of these suggested packages must be installed.  Various problems are theoretically expected for \eqn{p} close to zero (see Details), so the default parameter range when \eqn{p} is estimated is 1.1 -- 1.999.

A structured-dispersion model can be fitted for \eqn{\phi} but this is sensitive to numerical issues (see Details about \dQuote{inner estimation} of \eqn{\phi}). 

}
\usage{
tweedie(index, link = "log", numderiv=TRUE, ...)
}
\arguments{
  \item{index}{
Missing or numeric value for the index parameter \eqn{p}, subject to the constraint \eqn{1<p<2}. The \dQuote{outer-optimizing} functions, such as \code{fitme} (see \code{\link{inner-vs-outer}}) will estimate this parameter if called with \code{family=tweedie()} (no given \code{index} value). In other usages of this family, \code{index} must be specified.
}
  \item{link}{
Link function, specified either as the default \code{"log"} for a log link, or as the link exponent \eqn{q}. \eqn{q=0} is interpreted as the log link. There is no definite constraint on \eqn{q}, but using values that imply that \eqn{\eta} should be positive (such as 1, or non-integer values) may run into troubles. Values closes to 0 and 1 also lead to numerical errors. See Details if you want to estimate this parameter.  
}
  \item{numderiv}{
Boolean; whether to use numerical derivation in some computation. \code{numderiv=FALSE} (not recommended) will be ignored if the \pkg{tweedie} package is not installed.    
}
  \item{\ldots}{
Additional arguments passed to some member functions of the returned family object. Currently used only for the \code{method} and \code{method.args} of \code{\link[numDeriv]{grad}} if numerical derivatives are used (see \code{numderiv} argument).    
}
}
\details{
One can force \pkg{spaMM}'s \dQuote{outer-optimization} fitting functions to estimate the \code{link} exponent by using \code{tweedie(link="fit")}, and specyfing init, lower and upper values for the link exponent (using the name \code{Tw_link} for their values: see Examples). However this will fail if negative values of \eqn{\eta} are generated during the fitting process.

Estimates of the residual dispersion parameter \eqn{\phi} are by default constrained to be > 1e-05 because evaluations of the Tweedie density function are increasingly slow as \eqn{\phi} approaches 0. Other minimal \eqn{\phi} value can be set through the \code{lower} argument, as usual (or through \code{\link{control.HLfit}$min_phi} for inner estimation).

When \eqn{p} is close to 1, the Tweedie density is multimodal; further, when responses are recorded to a limited number of decimal places, the profile likelihood for \eqn{p} may be bimodal with a spurious maximum at \eqn{p = 1} (Dunn & Smyth 2005). 
From limited numerical experiments, for low \eqn{p} such as 1.05, computed likelihood profiles with respect to \eqn{\phi} also often exhibited several maxima. Overall, estimation for low \eqn{p} will be sensitive to many things.\cr

If the fitting function uses inner-estimation methods to fit \eqn{\phi} (see \code{\link{inner-vs-outer}}), it then requires evaluation of a derivative of the Tweedie density, which may be slow or not very accurate, depending on the method used (which can be controlled by the \code{numderiv} and optional \ldots arguments). So the fitting function avoids inner estimation for \eqn{\phi}, when possible. When it is used, \cr
\eqn{\bullet} examples were found where the inner-estimation iterative algorithm did not converge due to cyclic values of \eqn{\phi} estimates, for \eqn{p<1.5}. An experimental correction has been implemented to fix this convergence issue;\cr
\eqn{\bullet} the results are sensitive to the algorithm used to compute the derivative of the conditional likelihood (with the default algorithm appearing better), and maybe more so when both \eqn{p} and \eqn{\phi} are small.

}
\value{
A family object with structure similar to \code{stats::}\code{\link[stats]{family}} objects, but with additional member functions for usage with \pkg{spaMM} fitting functions.
}
\references{
Dunn, Peter K and Smyth, Gordon K (2005) Series evaluation of Tweedie exponential dispersion model densities. Stat Comput 15, 267–280. \doi{10.1007/s11222-005-4070-y}
  
Dunn, Peter K and Smyth, Gordon K (2018) Generalized Linear Models With Examples in R. Springer, NY.
  \doi{10.1007/978-1-4419-0118-7}
  
Dunn, P. K. (2026).  Tweedie: Evaluation of Tweedie exponential family models. R package.
 \doi{10.32614/CRAN.package.tweedie}
}
\examples{
if (spaMM.getOption("example_maxtime")>1.4) {
### toy data 
  set.seed(123)
  y <- rgamma(20,shape=rpois(20,lambda=rgamma(20,5)))
  x <- 1:20
  toyTw <- data.frame(x=x,y=y, grp=x)
  
### Fits
fitme(y~x+(1|grp),family=spaMM::tweedie(), data=toyTw) 

# tentative estimation of the link exponent:      
fitme(y~x+(1|grp),family=spaMM::tweedie(link="fit"), data=toyTw, 
      init=list(lambda=NA, phi=NA, Tw_link=0.5),
      lower=list(Tw_link=0.4), upper=list(Tw_link=0.6)) 
# (this shows that this _can_ work).

### Interactions with other packages
  tweedie::logLiktweedie(toyglm <- 
    glm(y~x,family=spaMM::tweedie(index=1.5,link= -0.5), data=toyTw))
  toyglm$fitted.values
  predict(toyglm, type="response")
  tweedie::logLiktweedie(spaMM_glm(y~x,
    family=spaMM::tweedie(index=1.5,link= -0.5), data=toyTw))
}
}
\keyword{models}
\keyword{regression}