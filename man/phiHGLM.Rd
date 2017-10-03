\name{phiHGLM}
\alias{phiHGLM}

\title{Fitting random effects in the residual dispersion model}
\description{
  \eqn{\phi} parameters are estimated by fitting a Gamma HGLM to response values computed by the parent fitting function (e.g., by \code{HLfit} in the Examples). The \code{fitme} function is used to perform this fit. The \code{resid.model} of the parent call is used to control the arguments of this \code{fitme} call. 
}
\usage{
# 'resid.model' argument of main fitting functions
}
\arguments{
  \code{resid.model} is \bold{either} a formula (without left-hand side) for the dispersion parameter \code{phi} of the residual error (a log link is assumed);\cr 
  \bold{or} a list, with possible elements:
    \item{formula}{model formula as in formula-only case, without left-hand side}
    \item{family}{The family is always Gamma. The default link is log. The identity link can be tried but may fail because only the log link ensures that the fitted \eqn{\phi} is positive.}
    % inverse cannot be tried
    \item{fixed}{fixed values of parameters. Same usage as documented in \code{\link{fitme}}}
    \item{control.dist}{A list of arguments that control the computation of the distance argument of the correlation functions. Same usage as documented in \code{\link{HLCor}}} 
    \item{rand.family}{A \code{family} object or a \code{list} of family objects describing the distribution of the random effect(s). Same usage as documented for \code{\link{HLfit}}} 
  \code{resid.model} with random effects is still experimental and complex combinations of arguments could give unexpected results. In particular, the functionality of \code{init.HLfit, lower, upper, control} has not been tested. 
  The list \code{should not}  contain the following elements:
    \item{init}{Currently ignored;}
    \item{method}{which is constrained to be identical to the method from the parent call;} 
    \item{control.HLfit, control.glm}{constrained to be identical to the same-named controls from the parent call;}
    \item{resid.model}{(constrained: no \code{resid.model} for a \code{resid.model});}
    \item{REMLformula}{(constrained to NULL);}
    \item{data}{identical to data from the parent call, which must therefore include all the variables required for the \code{resid.model};}
    \item{prior.weights}{constrained: no prior weights;}
    \item{verbose}{constrained: will display a progress line summarizing the results of the \code{resid.model} fit at each iteration of main loop of the parent \code{HLfit} call.}
}
\references{
Lee, Y., Nelder, J. A. and Pawitan, Y. (2006) Generalized linear models with random effects: unified analysis via
h-likelihood. Chapman & Hall: London.
}
\examples{
if (spaMM.getOption("example_maxtime")>4.9) {
 data("crack") # crack data, Lee et al. 2006 chapter 11 etc
 hlfit <- HLfit(y~crack0+(1|specimen),family=Gamma(log),
                data=crack, HLmethod="REML", 
                rand.family=inverse.Gamma(log), 
                resid.model=list(formula=~cycle+(1|specimen)) )
}
}