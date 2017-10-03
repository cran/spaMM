\name{spaMM-conventions}
\alias{spaMM-conventions}

\title{spaMM conventions and differences from related fitting procedures}

\description{
\bold{input arguments} are generally similar to those of \code{glm} and \code{(g)lmer}, in particular for the \code{spaMM::fitme}
function, with the exception of the \code{prior.weights} argument, which is simply \code{weights} in the other packages.
The name \code{prior.weights} seems more consistent, since e.g. \code{glm} returns its input \code{weights} as output \code{prior.weights}, while its output \code{weights} are instead the weights in the final iteration of an iteratively weighted least-square fit.

The \\bold{default likelihood target} for dispersion parameters is restricted likelihood (REML estimation) for \code{corrHLfit} and (marginal) likelihood (ML estimation) for \code{fitme}. 
Model fits may provide restricted likelihood values(\code{ReL}) even if restricted likelihood is is not used as an objective function at any step in the analysis.   

See \code{\link{good-practice}} for advice about the proper syntax of \code{formula}.

\bold{Computation times} depend on control parameters given by  \code{spaMM.getOption("spaMM_tol")} parameters (for iterative algorithms), and \code{spaMM.getOption("nloptr")} parameters for the default optimizer. Do not use \code{spaMM.options()} to control them globally, unless you know what you are doing. Rather control them locally by the \code{control.HLfit} argument to control \code{spaMM_tol}, and by the control arguments of \code{corrHLfit} and \code{fitme} to control \code{nloptr}. If \code{nloptr$Xtol_rel} is set above 5e-06, \code{fitme} will by default refit the fixed effects and dispersion parameters (but not other correlation parameters estimated by \code{nloptr}) by the iterative algorithm after \code{nloptr} convergence. Increasing \code{nloptr$Xtol_rel} value may therefore switches the bulk of computation time from the optimizer to the iterative algorithm, and may increase or decrease computation time depending on which algorithm is faster for a given input. Use \code{control$refit} if you wish to inhibit this, but note that by default it provides a rescue to a poor \code{nloptr} result due to a too large \code{Xtol_rel}.       
}
\references{
Chambers J.M. (2008) Software for data analysis: Programming with R. Springer-Verlag New York
}