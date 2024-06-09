\name{hatvalues.HLfit}
\alias{hatvalues}
\alias{hatvalues.HLfit}
\title{
Leverage extractor for HLfit objects
}
\description{
This gets \dQuote{leverages} or \dQuote{hat values} from an object. However, there is hidden complexity in what this may mean, so care must be used in selecting proper arguments for a given use (see Details). To get the full hat matrix, see \code{\link{get_matrix}(., which="hat_matrix")}.
}
\usage{
\method{hatvalues}{HLfit}(model, type = "projection", which = "resid", force=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{An object of class \code{HLfit}, as returned by the fitting functions in \code{spaMM}.}
  \item{type}{Character: \code{"projection"}, \code{"std"}, or more cryptic values not documented here. See Details.} 
  \item{which}{Character: \code{"resid"} for the traditional leverages of the observations, \code{"ranef"} for random-effect leverages, or \code{"both"} for both.}
  \item{force}{Boolean: to force recomputation of the leverages even if they are available in the object, for checking purposes.}
  \item{\dots}{For consistency with the generic.}
}
\details{
Leverages may have distinct meaning depending on context. The textbook version for linear models is that leverages \eqn{(q_i)} are the diagonal elements of a projection matrix (\dQuote{hat matrix}), and that they may be used to standardize (\dQuote{studentize}) residuals as follows. If the residual variance \eqn{\phi} is known, then the variance of each fitted residual \eqn{\hat{e}_i} is \eqn{\phi(1-q_i)}. Standardized residuals, all with variance 1, are then \eqn{\hat{e}_i/}\eqn{\sqrt{}}{\sqrt}\eqn{(\phi(1-q_i))}. This standardization of variance no longer holds exactly with estimated \eqn{\phi}, but if one uses here an unbiased (REML) estimator of \eqn{\phi}, the studentized residuals may still practically have a unit expected variance. 

When a simple linear model is fitted by ML, the variance of the fitted residuals is less than \eqn{\phi}, but \eqn{\hat{\phi}} is downward biased so that residuals standardized only by \eqn{\sqrt{}}{\sqrt}\eqn{(\phi)}, without any leverage correction, more closely have expected unit variance than if corrected by the previous leverages. The ML and REML computations can be seen as both using \dQuote{standardizing} leverages, defined so that they are zero in the ML case and are equal to the \dQuote{projection} leverages (the above ones, derived from a projection matrix) in the REML case.

These \dQuote{standardizing} leverages can themselves been seen as special cases of those that appear in expressions for derivatives, with respect to the dispersion parameters, of the log-determinant of the information matrices considered in the Laplace approximation for marginal or restricted likelihood (Lee et al. 2006). This provides a basis to generalize the concept of standardizing leverages for ML and REML in mixed-effect models. In particular, in an ML fit, one considers leverages \eqn{(q*_i)} that are no longer the diagonal elements of the projection matrix for the mixed model [and, as hinted above, for a simple linear model the ML \eqn{(q*_i)} are zero]. The generalized standardizing leverages may include corrections for non-Gaussian response, for non-Gaussian random effects, and for taking into account the variation of the GLM weights in the logdet(info.mat) derivatives. Which corrections are included depend on the precise method used to fit the model (e.g., EQL vs PQL vs REML). Standardizing leverages are also defined for the random effects.

These distinctions suggest breaking the usual synonymy between \dQuote{leverages} or \dQuote{hat values}: the term \dQuote{hat values} better stands for the diagonal elements of a projection matrix, while \dQuote{leverages} better stands for the standardizing values.   
\code{hatvalues(.,type="std")} returns the standardizing leverages. By contrast, \code{hatvalues(.,type="projection")} will always return hat values from the fitted projection matrix. Note that these values typically differ between ML and REML fit because the fitted projection matrix differs between them.
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
