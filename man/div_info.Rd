\name{div_info}
\alias{div_info}
\alias{diagnose_conv}
\title{
Information about numerical problems
}
\description{
This experimental function displays information about parameter values for which some numerical problems have occurred. Some warnings suggest its use.

Numerical problems may occur if the information matrix (for the augmented linear model used in the iteratively reweighted least-squares algorithm) is nearly singular. \pkg{spaMM} may try to check whether such singularity occurs when this algorithm has not converged. But itself may be slow so it is not performed systematically for large matrices. \code{spaMM.options(diagnose_conv=<integer>)} may be used to control the maximum size of matrices for which the check is performed.

When \dQuote{outer} generic optimization is performed, information is reported about the range of parameter values for which problems occurred, (see Value). The fit object \code{divinfo} element may also contain more informative tables of parameter values. This information is currently missing for \dQuote{inner}-optimized parameters.
}
\usage{
div_info(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of class \code{HLfit}, as returned by the fitting functions in \code{spaMM}.}
  \item{\dots}{Currently not used}
}
\value{Used mainly for the side effects (printed output) but returns invisibly either a single parameter vector (if a single numerical problem occurred) or a matrix of parameter ranges, or NULL if there is no problem to report.
}
\examples{
# Tragically ;-), no simple example of numerical problems % ____F I X M E____
# that can be diagnosed by div_info() is currently available.
}
