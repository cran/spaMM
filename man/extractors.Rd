\name{extractors}
\alias{model.frame.HLfit} % undocumented as long as it is not stable (this has limitations)
\alias{model.matrix.HLfit} % FIXME not used in spaMM code, and undocumented
\alias{extractors}
\alias{getDistMat}
\alias{logLik}
\alias{logLik.HLfit}
\alias{fitted}
\alias{fitted.HLfit}
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
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Functions to extract various components of a fit
}
\description{
\code{formula} extracts the model formula.
\code{family} extracts the response family.
\code{terms} extracts the formula, with attributes describing the \bold{fixed-effect} terms. 
\code{nobs} returns the length of the response vector.
\code{logLik} extracts the log-likelihood (exact or approximated).
\code{dev_resids} returns a vector of squared (unscaled) deviance residuals (the summands in McCullagh and Nelder 1989, p. 34). 
\code{deviance} returns the sum of squares of these (unscaled) deviance residuals, that is (consistently with \code{stats::deviance}) the unscaled deviance.
\code{fitted} extracts fitted values (see \code{\link{fitted.values}}).
\code{response} extracts the response (as a vector).
\code{fixef} extracts the fixed effects coefficients, \eqn{\beta}.
\code{ranef} extracts the predicted random effects, \bold{Lv} (default since version 1.12.0), or \bold{u} (see Details in \code{\link{HLfit}} for definitions), \code{print.ranef} controls their printing. 
\code{getDistMat} returns a distance matrix for a geostatistical (Matérn etc.) random effect.
}
\usage{
\method{formula}{HLfit}(x, which="hyper", ...)
\method{family}{HLfit}(object, ...)
\method{terms}{HLfit}(x, ...)
\method{nobs}{HLfit}(object, ...)
\method{logLik}{HLfit}(object, which, ...)
\method{fitted}{HLfit}(object, ...)
\method{fixef}{HLfit}(object, ...)
\method{ranef}{HLfit}(object, type = "correlated", ...)
\method{print}{ranef}(x, max.print = 40L, ...)
\method{deviance}{HLfit}(object, ...)
##
getDistMat(object, scaled=FALSE, which = 1L)
response(object,...)
dev_resids(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of class \code{HLfit}, as returned by the fitting functions in \code{spaMM}.}
  \item{type}{For \code{ranef}, use \code{type="correlated"} (default) to display the correlated random effects (\bold{Lv}), whether in a spatial model, or a random- coefficient model. Use \code{type="uncorrelated"} to pretty-print the elements of the \code{<object>$ranef} vector (\bold{u}).}
  \item{which}{For \code{logLik}, the name of the element of the \code{APHLs} list to return (see Details for any further possibility). The default depends on the fitting method. In particular, if it was REML or one of its variants, the function returns the log restricted likelihood (exact or approximated). For \code{getDistMat}, an integer, to select a random effect from several for which a distance matrix may be constructed. For \code{formula}, by default the model formula with non-expanded \code{multIMRF} random-effect terms is returned, while for \code{which=""} a formula with \code{multIMRF} terms expanded as \code{IMRF} terms is returned.}
  \item{scaled}{If \code{FALSE}, the function ignores the scale parameter \eqn{rho} and returns unscaled distance.}
  \item{x}{For \code{print.ranef}: the return value of \code{ranef.HLfit}.}
  \item{max.print}{Controls \code{options("max.print")} locally.}
  \item{\dots}{Other arguments that may be needed by some method.}
}
\value{
\code{formula} returns a \code{formula}, except a list of them from \code{fitmv()} output.

\code{terms} returns an object of class \code{c("terms", "formula")} which contains the \emph{terms} representation of a symbolic model.  See \code{\link[stats]{terms.object}} for its structure. \code{terms(<fitmv() result>)} returns a list of such terms.

Other return values are numeric (for \code{logLik}), vectors (most cases), matrices or dist objects (for \code{getDistMat}), or a family object (for \code{family}). \code{ranef} returns a list of vectors or matrices (the latter for random-coefficient terms). 
}
\details{
 With \code{which="LogL_Lap"}, \code{logLik()} returns a Laplace approximation of log-likelihood based on the observed Hessian, rather than the expected Hessian. This is implemented only for the case \code{family=Gamma(log)}, for demonstration purposes.
}
\seealso{
See \code{\link{residuals.HLfit}} to extract various residuals, \code{\link{residVar}} to extract residual variances or information about residual variance models, \code{\link{hatvalues}} to extract leverages, \code{\link{get_matrix}} to extract the model matrix and derived matrices, and \code{\link{vcov.HLfit}} to extract covariances matrices from a fit, \code{\link{get_RLRsim_args}} to extract arguments for (notably) tests of random effects in LMMs.
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
m1 <- HLfit(y ~ X1+X2+(1|batch), resid.model = ~ 1,
            data=wafers, method="ML")
fixef(m1)
ranef(m1)

data("blackcap")
fitobject <- fitme(migStatus ~ 1 + Matern(1|longitude+latitude),data=blackcap,
                       fixed=list(nu=4,rho=0.4,phi=0.05))
getDistMat(fitobject)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
