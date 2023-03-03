\name{phi-resid.model}
\alias{phi-resid.model}

\title{Residual dispersion model for gaussian and Gamma response}
\description{
A model can be specified for the residual-dispersion parameter \eqn{\phi} of gaussian and Gamma response families. The \code{resid.model} argument of all fitting functions is used to specify this model. This model may or may not include random effects. It is fitted by a specific method (Lee & Nelder 2006) involving estimation of its parameters by the fit of a Gamma-response model to response values computed by the parent fitting function (e.g., by \code{HLfit} in the Examples). The \code{fitme} function is used internally to perform this fit, irrespective of the parent fitting function. The \code{resid.model} argument of the parent call is used to control the arguments of this \code{fitme} call beyond the specification of the model. 

When the residual-dispersion model includes random effects, no single likelihood objective function appears to be maximized by the joint fit of mean-response and residual dispersion models. A procedure such as \code{numInfo} may then detect that the likelihood gradient does not vanish for all parameters. Indeed, this limitation is \dQuote{relatively obvious} in Lee & Nelder's original formulation since (as classical REML methods do) they used marginal likelihood and restricted likelihood concepts to fit different parameters of the joint model. But this limitation is also true in the case where marginal likelihood (actually, its Laplace approximation, although the issue could persist even if exact Gamma-GLMM likelihood were used) is used in the residual-dispersion fit.     
}
\usage{
# 'resid.model' argument of fitting functions (fitme(), HLfit(), etc)
}
\arguments{
  \code{resid.model} is \bold{either} a formula (without left-hand side) for the dispersion parameter \code{phi} of the residual error (a log link is assumed); \bold{or} a list of arguments similar to those of a standard fit. The following arguments may be useful:
  %
  \item{formula}{model formula as in formula-only case, without left-hand side.}
  \item{family}{The family is always Gamma. The default link is log. The identity link can be tried but may fail because only the log link ensures that the fitted \eqn{\phi} is positive.
        % inverse cannot be tried
        }
  \item{fixed}{fixed values of parameters of the residual dispersion model itself. Same usage as documented in \code{\link{fitme}}, except that it is better not to try to fix its \code{phi} (see Details).} 
  \item{etaFix}{To fix some of the fixed-effect coefficients, as in the mean response, and with the same format. Note that the same effect can usually be acheived by an offset in the \code{formula}.} 

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
    \item{init.HLfit}{if used, this argument may affect the fits. However, it is best ignored in practice: users would have hard time guessing good initial values, and they might have unwarranted effects.}
  }
  the \code{phi} of the Gamma family of the residual dispersion model is by default set to 1, in agreement with the theory underlying the estimation procedure for the residual model; it can be set to another value, and a \code{resid.model}'s \code{fixed=list(phi=NA)} will even force its estimation, but this is not warranted.

Fits with a mixed-effect residual-dispersion model involve repeated (\dQuote{fitme} fits of the latter model (themselves within the \dQuote{HLfit} calls nested within the main fit), which can be slow particularly when this model involve spatial effects. A specific element \code{phifit} of the \code{verbose} vector controls screen information about progress of such fits during the full-model fit: when set to 0 (or FALSE) there is no report. For higher values a one-line message is output at the end of each nested \dQuote{fitme} call, but it may be overwritten by the next one line message. So the ultimately visible output depends on control of overwriting. When \code{verbose["phifit"]} is set to 1 (or TRUE) each output overwrites the previous one so the ultimately visible output is from the last nested \dQuote{fitme} call; when it is set to \code{2}, theline of output of the final nested \dQuote{fitme} call remains visible for each \dQuote{HLfit} call; when set to 3, a line of output remains visible from each nested \dQuote{fitme} call within each \dQuote{HLfit} call.  
}
\value{
When such dispersion models are fitted, the resulting fits are embedded in the main fit object. The \code{\link{get_fittedPars}} extractor will by default )as controlled by its argument \code{phiPars}) include in its return value the \code{rdisPars} element, which is the list of parameters of the residual-dispersion fit, in the same format as a \code{get_fittedPars} value for the mean-response model (\code{rdisPars} may also include fits of other residual-dispersion models described in \code{\link{resid.model}}). The \code{phi} element of the \code{get_fittedPars} value will further contain the residual-dispersion fit itself, as a \code{"glm"} or, when it includes random effects, as a \code{"HLfit"} object.
}
\references{
Lee, Y. and Nelder, J.A. (2006), Double hierarchical generalized linear models (with discussion). Journal of the Royal Statistical Society: Series C (Applied Statistics), 55: 139-185. \doi{10.1111/j.1467-9876.2006.00538.x}

Lee, Y., Nelder, J. A. and Pawitan, Y. (2006) Generalized linear models with random effects: unified analysis via
h-likelihood. Chapman & Hall: London.
}
\examples{
 data("crack") # crack data, Lee et al. 2006 chapter 11 etc
 hlfit <- HLfit(y~crack0+(1|specimen), family=Gamma(log),
                data=crack, rand.family=inverse.Gamma(log), 
                resid.model=list(formula=~cycle+(1|specimen)) )
}