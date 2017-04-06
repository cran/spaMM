\name{extractors}
\alias{extractors}
\alias{getDistMat}
\alias{logLik}
\alias{logLik.HLfit}
\alias{fitted}
\alias{fitted.HLfit}
\alias{residuals}
\alias{residuals.HLfit}
\alias{fixef}
\alias{fixef.HLfit}
\alias{ranef}
\alias{ranef.HLfit}
\alias{vcov}
\alias{vcov.HLfit}
\alias{deviance}
\alias{deviance.HLfit}
\alias{dev_resids}
\alias{get_any_IC}
\alias{AIC}
\alias{AIC.HLfit}
\alias{Corr}
\alias{get_ZALMatrix}
\alias{get_RLRTSim_args}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Functions to extract various components of a fit
}
\description{
\code{logLik} extracts the log-likelihood (exact or approximated).
\code{get_any_IC} and \code{AIC.HLfit} compute model selection/information criteria such as AIC. See \code{\link{AIC}} for more information about these criteria.
\code{dev_resids} returns a vector of deviance residuals. 
\code{deviance} returns the sum of these deviance residuals.
\code{fitted} extracts fitted values (see \code{\link{fitted.values}}).
\code{residuals} extracts residuals of the fit.
\code{fixef} extracts the fixed effects coefficients, \eqn{\beta}.
\code{ranef} extracts the predicted random effects, \bold{Lv} (default since version 1.12.0), or \bold{u} (see Details in \code{\link{HLfit}} for definitions).
\code{vcov} returns the variance-covariance matrix of the fixed-effects coefficients.
\code{Corr} returns a correlation matrix of random effects.
\code{getDistMat} extracts a distance matrix for a Matérn correlation model.
\code{get_ZALMatrix} extracts the design matrix for the random effects \eqn{v}.
\code{get_RLRTSim_args} extracts a list of arguments suitable for calls to \code{LRTSim::RLRTSim()}
}
\usage{
\method{logLik}{HLfit}(object,which,...)
\method{fitted}{HLfit}(object,...)
\method{fixef}{HLfit}(object,...)
\method{ranef}{HLfit}(object, type="correlated", ...)
\method{vcov}{HLfit}(object,...)
\method{deviance}{HLfit}(object,...)
getDistMat(object,scaled=FALSE)
Corr(object,...)
dev_resids(object,...)
get_any_IC(object, ..., verbose=interactive())
get_RLRTSim_args(object,...)
get_ZALMatrix(object,as_matrix)
\method{AIC}{HLfit}(object, ..., k, verbose=interactive())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{The return object of an HLfit or similar function.}
  \item{type}{Use \code{type="correlated"} (default) to display the correlated random effects (\bold{Lv}), whether in a spatial model, or a random- coefficient model. Use \code{type="uncorrelated"} to pretty-print the elements of the \code{<object>$ranef} vector (\bold{u}).}
  \item{which}{Which element of the \code{APHLs} list to return. The default depends on the fitting method.In particular, if it was REML or one of its variants, the function returns the log restricted likelihood (exact or approximated).}
  \item{scaled}{If \code{FALSE}, the function ignores the scale parameter \eqn{rho} and returns unscaled distance.}
  \item{as_matrix}{Boolean: whether to always return a \code{matrix}; otherwise as \code{Matrix} may be returned.}
  \item{verbose}{ Whether to print the model selection criteria or not. }
  \item{k}{For \code{AIC}, unused by \code{HLfit} method, but included to conform to the generic.}
  \item{\dots}{Other arguments that may be needed by some method.}
}
\value{
Return values are numeric (for \code{logLik}), vectors (most cases), matrices (for \code{vcov}), matrices or dist objects (for \code{getDistMat}). \code{ranef} returns a vector with attributes, which inherits from class \code{ranef} which has its own (undocumented) \code{print} method. 

\code{Corr} returns a list, for the different random effect terms, of unconditional correlation matrix of the random effects \dQuote{\bold{v}} (see Details of \code{\link{HLfit}} for definitions).

\code{get_any_IC} computes, optionally prints, and returns invisibly the following quantities. The \bold{Effective degrees of freedom} for the random effects (approximately) characterizes the expectation of a goodness of fit statistic discussed by Lee and Nelder (2001), which gave a general formula for it in HGLMs. The \bold{conditional AIC} (Vaida and Blanchard 2005) is notable in involving the conditional likelihood and the effective degrees of freedom. Lee et al. (2006) and Ha et al (2007) defined a corrected AIC [i.e., AIC(D*) in their eq. 7] that is here interpreted as the conditional AIC. The conditional AIC returned by \code{HLfit} includes both this effective df, the df for estimated fixed effects, and the df for estimated parameters of the variance of random effects. 
% FR->FR If \code{HLfit} was called through \code{corrHLfit}, a correction for estimated correlation parameters may be added.  
Also returned are the \bold{marginal AIC} (Akaike's classical AIC), and a focussed AIC for dispersion parameters (\bold{dispersion AIC}) discussed by Ha et al (2007; eq.10). This diversity of criteria should encourage users to think twice before applying model selection automatically, which is no better although more fashionable than misuses of simple null hypothesis testing. Also, alternative procedures for model choice can be considered (e.g. Cox and Donnelly, 2011, p. 130-131).

\code{get_RLRTSim_args} extracts a list of arguments suitable for a call to \code{LRTSim::RLRTSim()} for a small-sample test of the presence of a random effect by an efficient simulation procedure. The test can be run by\cr \code{do.call("RLRTSim",<get_RLRTSim_args return value>)}.  
}
\references{
Cox, D. R. and Donnelly C. A. (2011) Principles of Applied Statistics. Cambridge Univ. Press.

Ha, I. D., Lee, Y. and MacKenzie, G. (2007) Model selection for multi-component frailty models. Statistics in Medicine 26: 4790-4807.

Lee, Y., Nelder, J. A. (2001)  Hierarchical generalised linear models: A
synthesis of generalised linear models, random-effect models and structured
dispersions. Biometrika 88, 987-1006.

Lee, Y., Nelder, J. A. and Pawitan, Y. (2006) Generalized linear models with random effects: unified analysis via
h-likelihood. Chapman & Hall: London.

Vaida, F., and Blanchard, S. (2005) Conditional Akaike information for mixed-effects models. Biometrika 92, 351-370.
}
\examples{
data(wafers)
m1 <- HLfit(y ~X1+X2+(1|batch),
          resid.model = ~ 1 ,data=wafers,HLmethod="ML")
get_any_IC(m1)          
fixef(m1)
vcov(m1)
ranef(m1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
