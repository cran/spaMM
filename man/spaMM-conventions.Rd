\name{spaMM-conventions}
\alias{spaMM-conventions}

\title{Difference between spaMM and common R practices}

\description{
\code{spaMM} differs in the following ways from practices one may be used to:

\bold{input arguments} are generally similar to those of \code{glm} and \code{(g)lmer}, in particular for the \code{spaMM::fitme}
function, with the exception of the \code{prior.weights} argument, which is simply \code{weights} in the other packages.
\code{prior.weights} is in a sense more consistent, since e.g. \code{glm} returns its input \code{weights} a output \code{prior.weights}, while its output \code{weights} are the distinct GLM weights.

The default \bold{likelihood target} for dispersion parameters is restricted likelihood (REML estimation) for \code{corrHLfit} and (marginal) likelihood (ML estimation) for \code{fitme}. 

Most basic fitting functions in R will seek variables in the global environment if they are not in the \code{data}. This easily leads to errors (see example in the discussion of \code{\link{update.HLfit}}) and \code{spaMM} tries to restrain this.

\bold{Computation times} depend on \code{control.HLfit$conv.threshold}, which is by default 1e-5. Increasing it to 1e-4 had no notable effect (08/2016) on the tests in the \code{tests} directory of the package, except reducing computation time. \code{spaMM} is therefore slower but safer than it could be.     
}
