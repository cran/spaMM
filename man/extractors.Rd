\name{extractors}
\alias{model.frame.HLfit} % undocumented as long as it is not stable (and not yet used)
\alias{model.matrix.HLfit} % FIXME not used in spaMM code, and undocumented
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
\alias{formula}
\alias{formula.HLfit}
\alias{family}
\alias{family.HLfit}
\alias{terms}
\alias{terms.HLfit}
\alias{nobs}
\alias{nobs.HLfit}
\alias{ranef}
\alias{ranef.HLfit}
\alias{print.ranef}
\alias{deviance}
\alias{deviance.HLfit}
\alias{response}
\alias{dev_resids}
\alias{get_RLRTSim_args}
\alias{get_boot_response}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Functions to extract various components of a fit
}
\description{
\code{formula} extracts the model formula.
\code{family} extracts the response family.
\code{terms} extracts the \bold{fixed-effect} terms.
\code{nobs} returns the length of the response vector.
\code{logLik} extracts the log-likelihood (exact or approximated).
\code{dev_resids} returns a vector of squared (unscaled) deviance residuals (the summands in McCullagh and Nelder 1989, p. 34). 
\code{deviance} returns the sum of squares of these (unscaled) deviance residuals, that is (consistently with \code{stats::deviance}) the unscaled deviance.
\code{fitted} extracts fitted values (see \code{\link{fitted.values}}).
\code{residuals} extracts residuals of the fit.
\code{response} extracts the response (as a vector).
\code{fixef} extracts the fixed effects coefficients, \eqn{\beta}.
\code{ranef} extracts the predicted random effects, \bold{Lv} (default since version 1.12.0), or \bold{u} (see Details in \code{\link{HLfit}} for definitions), \code{print.ranef} controls their printing. 
\code{getDistMat} returns a distance matrix for a Matérn correlation model.
\code{get_RLRTSim_args} returns a list of arguments suitable for calls to \code{RLRsim::RLRTSim()}
\code{get_boot_response} returns the response vector of a given bootstrap replicate in an obsolete use of \code{fixedLRT} (see Details). 
}
\usage{
\method{formula}{HLfit}(x, ...)
\method{family}{HLfit}(object, ...)
\method{terms}{HLfit}(x, ...)
\method{nobs}{HLfit}(object, ...)
\method{logLik}{HLfit}(object,which,...)
\method{fitted}{HLfit}(object,...)
\method{fixef}{HLfit}(object,...)
\method{ranef}{HLfit}(object, type="correlated", ...)
\method{print}{ranef}(x, max.print=40L, ...)
\method{deviance}{HLfit}(object,...)
\method{residuals}{HLfit}(object, type = c("deviance", "pearson", "response"), ...)
getDistMat(object, scaled=FALSE, which=1L)
response(object,...)
dev_resids(object,...)
get_RLRTSim_args(object,...)
get_boot_response(object, replicate)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A object of class \code{HLfit}, as returned by the fitting functions in \code{spaMM}.}
  \item{type}{For \code{ranef}, use \code{type="correlated"} (default) to display the correlated random effects (\bold{Lv}), whether in a spatial model, or a random- coefficient model. Use \code{type="uncorrelated"} to pretty-print the elements of the \code{<object>$ranef} vector (\bold{u}). For \code{residuals}, the type of residuals which should be returned. The alternatives are: "deviance" (default), "pearson", and "response".}
  \item{which}{For \code{logLik}, the name of the element of the \code{APHLs} list to return (see Details for any further possibility). The default depends on the fitting method. In particular, if it was REML or one of its variants, the function returns the log restricted likelihood (exact or approximated). For \code{getDistMat}, an integer, to select a random effect from several for which a distance matrix may be constructed.}
  \item{scaled}{If \code{FALSE}, the function ignores the scale parameter \eqn{rho} and returns unscaled distance.}
  \item{x}{For \code{print.ranef}: the return value of \code{ranef.HLfit}.}
  \item{max.print}{Controls \code{options("max.print")} locally.}
  \item{replicate}{integer; to get the response vector of the given bootstrap replicate.}
  \item{\dots}{Other arguments that may be needed by some method.}
}
\value{
Return values are numeric (for \code{logLik}), vectors (most cases), matrices or dist objects (for \code{getDistMat}), or a family object (for \code{family}). \code{ranef} returns a list of vectors or matrices (the latter for random-coefficient terms). \code{terms} returns an object of class \code{c("terms", "formula")} which contains the \emph{terms} representation of a symbolic model.  See \code{\link[stats]{terms.object}} for its structure.

\code{get_boot_response} aims to reproduce samples produced by fixedLRT() with argument \code{boot_fn=".eval_boot_replicates"}. However, the default \code{boot_fn} is now \code{spaMM_boot} wherein samples are simulated by \code{simulate(object,nsim = nsim,verbose=FALSE,resp_testfn=resp_testfn)} and, given the same state of the random number generator, running this should be enough to reproduce the samples. \code{get_boot_response} has no \code{resp_testfn} argument so it cannot always reproduce the samples simulated by \code{spaMM_boot}. It may soon be deprecated, as \code{boot_fn=".eval_boot_replicates"}.

\code{get_RLRTSim_args} extracts a list of arguments suitable for a call to \code{RLRsim::RLRTSim()} for a small-sample test of the presence of a random effect by an efficient simulation procedure. The test can be run by\cr \code{do.call("RLRTSim",<get_RLRTSim_args return value>)}.  
}
\details{
See \code{\link{residuals.glm}} for more information about the types of residuals.

With \code{which="LogL_Lap"}, \code{logLik()} returns a Laplace approximation of log-likelihood based on the observed Hessian, rather than the expected Hessian. This is implemented only for the case \code{family=Gamma(log)}, for demonstration purposes.
}
\seealso{
See \code{\link{get_matrix}} \code{\link{vcov.HLfit}} to extract covariances matrices from a fit. 
}
\references{
McCullagh, P. and Nelder J. A. (1989) Generalized linear models. Second ed. Chapman & Hall: London.

Lee, Y., Nelder, J. A. (2001)  Hierarchical generalised linear models: A
synthesis of generalised linear models, random-effect models and structured
dispersions. Biometrika 88, 987-1006.

Lee, Y., Nelder, J. A. and Pawitan, Y. (2006) Generalized linear models with random effects: unified analysis via
h-likelihood. Chapman & Hall: London.
}
\examples{
data("wafers")
m1 <- HLfit(y ~X1+X2+(1|batch),
          resid.model = ~ 1 ,data=wafers,HLmethod="ML")
fixef(m1)
ranef(m1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
