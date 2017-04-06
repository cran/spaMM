\name{spaMM-conventions}
\alias{spaMM-conventions}

\title{spaMM conventions and differences from related fitting procedures}

\description{
\bold{input arguments} are generally similar to those of \code{glm} and \code{(g)lmer}, in particular for the \code{spaMM::fitme}
function, with the exception of the \code{prior.weights} argument, which is simply \code{weights} in the other packages.
The name \code{prior.weights} seems more consistent, since e.g. \code{glm} returns its input \code{weights} as output \code{prior.weights}, while its output \code{weights} are instead the so-called GLM weights.

The \\bold{default likelihood target} for dispersion parameters is restricted likelihood (REML estimation) for \code{corrHLfit} and (marginal) likelihood (ML estimation) for \code{fitme}. 
Model fits may provide restricted likelihood values(\code{ReL}) even if restricted likelihood is is not used as an objective function at any step in the analysis.   

See \code{\link{good-practice}} for advice about the proper syntax of \code{formula}.

\bold{Computation times} depend on \code{control.HLfit$conv.threshold}, which is by default 1e-5. Increasing it to 1e-4 had no notable effect (08/2016) on the tests in the \code{tests} directory of the package, except reducing computation time. \code{spaMM} is therefore slower but safer than it could be.     
}
\references{
Chambers J.M. (2008) Software for data analysis: Programming with R. Springer-Verlag New York
}