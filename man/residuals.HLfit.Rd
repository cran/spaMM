\name{residuals.HLfit}
\alias{residuals}
\alias{residuals.HLfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract model residuals
}
\description{
Extracts several types of residuals from an object of class \code{HLfit}. Note that the default type (\code{"deviance"}) of returned residuals differs from the default of equivalent functions in base R.  
}
\usage{
\method{residuals}{HLfit}(object, 
  type = c("deviance", "pearson", "response", "working", 
           "RQR", "std_RQR", "std_dev_res", "std_dev_rt"), 
  force=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of class \code{HLfit}, as returned by the fitting functions in \code{spaMM}.}
  \item{type}{
  The type of residuals which should be returned. See Details for additional information. 
  }
  \item{force}{Boolean: to force recomputation of the \code{"std_dev_res"} residuals even if they are available in the object, for checking purposes.}
  \item{\dots}{For consistency with the generic.}
}
\details{
The first four types \code{"deviance"} (default), \code{"pearson"}, \code{"response"} are \code{"working"} are, for GLM families, the same that are returned by \code{\link[stats]{residuals.glm}}. \code{"working"} residuals may be returned only for fixed-effect models. 

The default type is \code{"deviance"} as for \code{\link[stats]{residuals.glm}}, but since this may be confusing to some users, a gentle one-time warning will be issued when the \code{object} is a GLMM.

GLMs may be formulated as exponential dispersion models, i.e., models where the probability or probability density of each observation \eqn{y} can be written in the form
\eqn{f(y)=b(y, \phi) \exp[-d(y, \mu)/(2\phi)]} where \eqn{d(y, \mu)} is the \emph{unit deviance} for the observation, independent of the dispersion parameter \eqn{\phi} (e.g., Dunn & Smyth, 2018). The sum over observations, \eqn{\sum_i d(y_i,\mu_i)} is called the \emph{deviance}, and \eqn{\sum_i d(y_i,\mu_i)/\phi} is called the \emph{scaled deviance}, by McCullagh & Nelder (1989, pp. 33-34). In the \pkg{spaMM} documentation the former deviance is called the \emph{unscaled} deviance to emphasize this distinction. 

Prior weights \eqn{w} may be included in \dQuote{unscaled} quantities even though \eqn{\phi} is not.
This is consistent with what base R does: unscaled quantities returned by \code{stats::residuals(<glm>)} and \code{stats::deviance(<glm>)} will be different for equivalent fits with different prior weights (as produced by \code{glm(., family=gaussian, weights=rep(2,nrow<data>))} versus the \code{glm} call without weights). 

\code{residuals(., type="deviance")} returns the vector of deviance residuals for the different observations (as \code{\link[stats]{residuals.glm}}) does. The deviance residual for an observation is the signed square root of the (prior-weighted) unit deviance for this observation, sign\eqn{(y-\mu)\sqrt{w d(y,\mu)}}. 

The extractor \code{\link{dev_resids}} is equivalent to the \code{family} object's \code{dev.resids} function. Hence, it does not really returns residuals but rather their squares, i.e. the (prior-weighted) unit deviances, because this is already so for families defined in base R. 

The \code{\link{deviance}} extractor returns the sum of these prior-weighted unit deviances, i.e., the unscaled, prior-weighted deviance. This is consistent with what \code{stats::deviance(<glm>)} does.

\code{residuals(., type="std_dev_rt")} returns the (prior-weighted) \dQuote{standardized deviance residuals} as defined by Lee et al. (2006, p.52), i.e. the deviance residuals divided by \eqn{\sqrt{\phi (1-q)}}, where \eqn{\phi} is the dispersion parameter of the distribution family (a vector of values, for heteroscedastic cases), and \eqn{q} is a vector of leverages given by \code{hatvalues(., type="std")} (see \code{\link{hatvalues}} for details about these specific standardizing leverages). They are thus scaled by \eqn{\phi}, and signed.

\code{residuals(.type="std_dev_res")} returns the signed \bold{squares} of such standardized residuals.

  Some definitions must be extended for non-GLM response families. as in the latter case, there may not be a concept of unit deviance, and there is no conceptual distinction between scaled and unscaled residuals since the residual dispersion parameter (if any) is not generally a scale factor \eqn{\phi}. This amounts to the convention that  \eqn{\phi=1}. The deviance residuals for non-GLMs are as defined in Details of \code{\link{llm.fit}}. They are analogous to the scaled ones for GLMs as they depend on any dispersion parameter. \code{"std_dev_res"} and \code{"std_dev_rt"} residuals are defined from them as detailed above for GLM response families, and from the additional convention that \eqn{\phi=1}. Pearson residuals and response residuals are defined as in \code{stats:::residuals.glm}.  
The \code{"working"} residuals are defined for each observation as \eqn{- [d \log(}clik\eqn{)/d \eta]/[d^2 \log(}clik\eqn{)/d \eta^2]} where \sQuote{clik} is the conditional likelihood.

\code{residuals(.type="RQR")} returns the randomized quantile residuals (Dunn & Smyth, 1996), used by \code{\link{gof}}. \code{type="std_RQR"} returns the RQR residuals divided by \eqn{\sqrt{1-q}}, following the same logic as for standardized deviance residuals.

}
\value{ A vector of residuals}
\references{
Dunn, K. P., and Smyth, G. K. (1996). Randomized quantile residuals. Journal of Computational and Graphical Statistics 5, 1-10.  

Dunn, P. K., & Smyth, G. K. (2018). Generalized Linear Models With Examples in R. Springer.

Lee, Y., Nelder, J. A. and  Pawitan, Y. (2006). Generalized linear models with random effects: unified analysis via h-likelihood. Chapman & Hall: London.

McCullagh, P. and Nelder J. A. (1989) Generalized linear models. Second ed. Chapman & Hall: London.
}
\examples{
data("wafers")
fit <- fitme(y ~X1+(1|batch) ,data=wafers, init=list(phi=NaN))  # : this 'init' 
#                 implies that standardized deviance residuals are saved in the 
#                 fit result, allowing the following comparison: 

r1 <- residuals(fit, type="std_dev_res") # gets stored value
r2 <- residuals(fit, type="std_dev_res", force=TRUE) # forced recomputation
if (diff(range(r1-r2))>1e-14) stop()

##### 
\dontrun{
glmfit <- glm(I(y/1000)~X1, family=gaussian(), data=wafers)
deviance(glmfit) #           3...  (a)
sum(residuals(glmfit)^2) #   3...  (b) 

# Same model, with different parametrization of residual variance 
glmfit2 <- glm(I(y/1000)~X1, family=gaussian(), data=wafers, weights=rep(2,198))
deviance(glmfit2) #          6...  (c)  
sum(residuals(glmfit2)^2) #  6...  (d)

# Same comparison but for HLfit objects:
spfit <- fitme(I(y/1000)~X1, family=gaussian(), data=wafers)
deviance(spfit) #            3...  (e)
sum(residuals(spfit)^2) #    3...  (f) ~ sum(abs(residuals(.,"std_dev_res")))*phi
sum(dev_resids(spfit)) #     3...                         
# Gaussian case: "RQR" residuals ~ ("deviance" residuals)/sqrt(phi)

spfit2 <- fitme(I(y/1000)~X1, family=gaussian(), data=wafers, prior.weights=rep(2,198))
deviance(spfit2) #           6...  (g) ~ (c,d) # post v4.2.0
sum(residuals(spfit2)^2) #   6...  (h) ~ (c,d) 
                                       ~ sum(abs(residuals(.,"std_dev_res")))*phi/pw
sum(dev_resids(spfit2)) #    3...  
# "RQR" residuals still ~ ("deviance" residuals)/sqrt(phi)

# Unscaled residuals should not depend on arbitrarily fixed residual variance:
spfit3 <- fitme(I(y/1000)~X1, family=gaussian(), data=wafers, fixed=list(phi=2),
                prior.weights=rep(2,198))
deviance(spfit3) #           6...  (i) ~ (g)
sum(residuals(spfit3)^2) #   6...  (k) ~ (h)
sum(dev_resids(spfit3)) #    3...                         

% _F I X M E__ need exemple showing lack of effect of pw for non-GLM family?
}
}
