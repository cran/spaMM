\name{phiHGLM}
\alias{phiHGLM}
\alias{resid.model}

\title{Fitting random effects in the residual dispersion model}
\description{
  \eqn{\phi} parameters are estimated by fitting a Gamma HGLM to response values computed by the parent fitting function (e.g., by \code{HLfit} in the Examples). The \code{fitme} function is used to perform this fit. The \code{resid.model} of the parent call is used to control the arguments of this \code{fitme} call. 
}
\usage{
# 'resid.model' argument of main fitting functions
}
\arguments{
  \code{resid.model} is \bold{either} a formula (without left-hand side) for the dispersion parameter \code{phi} of the residual error (a log link is assumed); \bold{or} a list of arguments similar to those of a standard fit. The following arguments may be useful:
  %
    \item{formula}{model formula as in formula-only case, without left-hand side. Random effects can be included and this appears to work well in simple cases (block effects, or geostatistical models) but has not been tested, or hardly so, for other cases. }
    \item{family}{The family is always Gamma. The default link is log. The identity link can be tried but may fail because only the log link ensures that the fitted \eqn{\phi} is positive.}
    % inverse cannot be tried
    \item{fixed}{fixed values of parameters. Same usage as documented in \code{\link{fitme}}}
    \item{control.dist}{A list of arguments that control the computation of the distance argument of the correlation functions. Same usage as documented in \code{\link{HLCor}}} 
    \item{rand.family}{A \code{family} object or a \code{list} of family objects describing the distribution of the random effect(s). Same usage as documented for \code{\link{HLfit}}} 
    \item{init, lower, upper, control}{with same usage as documented in \code{\link{fitme}}, may be at least partly heeded.} 
    Other arguments should be ignored (see Details).
}
\details{
  The following elements in \code{resid.model} should be ignored:
  \describe{
    \item{method}{which is constrained to be identical to the method from the parent call;} 
    \item{control.HLfit, control.glm}{constrained to be identical to the same-named controls from the parent call;}
    \item{resid.model}{constrained: no \code{resid.model} for a \code{resid.model};}
    \item{REMLformula}{constrained to NULL;}
    \item{data}{The data of the parent call are used, so they must include all the variables required for the \code{resid.model};}
    \item{prior.weights}{constrained: no prior weights;}
    \item{verbose}{constrained: will display a progress line summarizing the results of the \code{resid.model} fit at each iteration of main loop of the parent call.}
    \item{init.HLfit}{Ignored. Users would have hard time guessing good initial values, which would be have to be ignored in most contexts anyway.}
  }
}
\references{
Lee, Y., Nelder, J. A. and Pawitan, Y. (2006) Generalized linear models with random effects: unified analysis via
h-likelihood. Chapman & Hall: London.
}
\examples{
 data("crack") # crack data, Lee et al. 2006 chapter 11 etc
 hlfit <- HLfit(y~crack0+(1|specimen), family=Gamma(log),
                data=crack, rand.family=inverse.Gamma(log), 
                resid.model=list(formula=~cycle+(1|specimen)) )
}