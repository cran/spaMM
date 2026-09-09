\name{hatvalues.HLfit}
\alias{hatvalues}
\alias{hatvalues.HLfit}
\title{
Hat values and leverages extractor
}
\description{
\code{hatvalues.HLfit} is a method for the generic function \code{\link[stats]{hatvalues}}:\cr
\eqn{\bullet} \code{hatvalues(.,type="projection")} will return \dQuote{hat values}, defined as diagonal elements of the projection matrix of the augmented linear model;\cr
\eqn{\bullet} \code{hatvalues(.,type="std")} will return generalized \dQuote{leverages}, which are standardizing valeus that can be used for estimation of dispersion parameters (see Details for the differences).    

To get the full hat matrix, see \code{\link{get_matrix}(., which="hat_matrix")}.
}
\usage{
\method{hatvalues}{HLfit}(model, type = "projection", which = "resid", force=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{An object of class \code{HLfit}, as returned by the fitting functions in \code{spaMM}.}
  \item{type}{Character: \code{"projection"}, \code{"std"}, or more cryptic values not documented.} 
  \item{which}{Character: \code{"resid"} for the traditional leverages of the observations, \code{"ranef"} for random-effect leverages, or \code{"both"} for both.}
  \item{force}{Boolean: to force recomputation of the leverages even if they are available in the object, for checking purposes.}
  \item{\dots}{For consistency with the generic.}
}
\details{
Leverages may have distinct meanings depending on context. The textbook version for linear models is that leverages \eqn{(q_i)} are the diagonal elements of a projection matrix (\dQuote{hat matrix}), and that they may be used to standardize (\dQuote{studentize}) residuals as follows. If the residual variance \eqn{\phi} is known, then the variance of each fitted residual \eqn{\hat{e}_i} is \eqn{\phi(1-q_i)}. Standardized residuals, all with variance 1, are then \eqn{\hat{e}_i/}\eqn{\sqrt{}}{\sqrt}\eqn{(\phi(1-q_i))}. This standardization of variance no longer holds exactly with estimated \eqn{\phi}, but if one uses here an unbiased (REML) estimator of \eqn{\phi}, the studentized residuals may still practically have a unit expected variance. 

In mixed models, a generalized definition of the projection matrix is available, in terms of the matrix of the augmented linear model [its value \eqn{\mathbf{P}_\mathrm{a}} can be retrieved by \code{\link{get_matrix}(., which="hat_matrix")}]. The \dQuote{hat values}, defined as the diagonal elements of \eqn{\mathbf{P}_\mathrm{a}}, are also defined for the random effects. However, in general the standardizing \dQuote{leverages} \eqn{\mathbf{q}^*} that can again be used to obtain estimates of \eqn{\phi} and \eqn{\lambda}, by standardizing residuals (for \eqn{\phi}), or by standardizing predicted random effects (for \eqn{\lambda}), are no longer simply the \dQuote{hat values}. 

For illustration, the following simple differences hold between \code{"projection"} hat values and \code{"std"} values:\cr
\eqn{\bullet} \code{"std"} values are defined so that they are 0 for ML fits of simple linear model.\cr
\eqn{\bullet} In LMMs fitted by REML, \code{"std"} values  are the same as \code{"projection"} values. For ML fits, these two forms of leverages differ, because \code{"std"} values are diagonal values of another matrix than \eqn{\mathbf{P}_\mathrm{a}}. 
Note that the \code{"projection"} hat values also differ between ML and REML fits because the fitted \eqn{\mathbf{P}_\mathrm{a}} itself depends on \eqn{\hat{\phi}} and  \eqn{\hat{\lambda}}, which differ between ML and REML fits. \cr
\eqn{\bullet} In GLMMs (or HGLMs), the \code{"std"} values may differ from the \code{"projection"} ones even for REML, because they may include corrections for non-Gaussian response family, for non-Gaussian random effects, and for taking into account the variation of the GLM weights in the logdet(info.mat) derivatives. Which corrections are included depend on the precise method used to fit the model (e.g., EQL vs PQL vs REML). 

}
\value{
A list with separate components \code{resid} (leverages of the observations) and \code{ranef} if \code{which="both"}, and a vector otherwise.
}
\references{
Lee, Y., Nelder, J. A. and Pawitan, Y. (2006) Generalized linear models with random effects: unified analysis via
h-likelihood. Chapman & Hall: London.
}
\examples{
if (spaMM.getOption("example_maxtime")>0.8) {
data("Orthodont",package = "nlme")
rnge <- (107:108)

# all different:
#
hatvalues(rlfit <- fitme(distance ~ age+(age|Subject), 
                         data = Orthodont, method="REML"))[rnge]
hatvalues(mlfit <- fitme(distance ~ age+(age|Subject), 
                         data = Orthodont))[rnge] 
hatvalues(mlfit,type="std")[rnge]
}
}
