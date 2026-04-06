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

\code{"deviance"} residuals are the signed square root of those returned by \code{\link{dev_resids}} when there are no prior weights. Note that there is a similar discrepancy in base R, as e.g. \code{gaussian()$dev.resids} returns the square of the \code{"deviance"} residuals. A similar discrepancy also affects the \code{"std_dev_res"} type.  Lee et al. (2006, p.52) defined the standardized deviance residuals as the deviance residuals divided by \eqn{\sqrt{\phi (1-q)}}, where \eqn{\phi} is the dispersion parameter of the distribution family (a vector of values, for heteroscedastic cases), and \eqn{q} is a vector of leverages given by \code{hatvalues(., type="std")} (see \code{\link{hatvalues}} for details about these specific standardizing leverages). \code{residuals(.type="std_dev_res")} returns the signed \bold{squares} of such standardized residuals. Use \code{residuals(.type="std_dev_rt")} to obtain the signed root values. 

The \code{"RQR"} type specifies the randomized quantile residuals (Dunn & Smyth, 1996), used by \code{\link{gof}}. \code{"std_RQR"} specifies the RQR residuals divided by \eqn{\sqrt{1-q}}, following the same logic as for standardized deviance residuals.

In the presence of prior weights, what the standard extractors do is often a matter of confusion and \pkg{spaMM} has not always been consistent with them. For a gaussian-response GLM (see Examples) \code{stats::deviance.lm} calls \code{\link{weighted.residuals}()} which returns \emph{unscaled} deviance residuals weighted by prior weights. Unscaled deviance residuals are defined in McCullagh and Nelder 1989, p. 34 and depend on the response values and fitted values but not on the canonical \eqn{\phi} parameter, and prior weights are not considered. \code{weighted.residuals()} still ignores \eqn{\phi} but accounts for prior weights. This means that different \code{residuals(<glm>)} and \code{deviance(<glm>)} will be returned for equivalent fits with different parametrizations of the residual variance (as produced by \code{glm(., family=gaussian, weights=rep(2,nrow<data>))} versus the \code{glm} call without weights). \code{residuals(<HLfit object>,"deviance")} and \code{deviance(<HLfit object>)} are consistent with this behavior. By contrast, \code{dev_resids(<HLfit object>)} always return the unscaled deviance residuals by default.   
Prior weights have an effect on the returned standardized deviance residuals in the same way as it has an effect on the deviance residuals.
  
  Some definitions must be extended for non-GLM response families. In the latter case, the deviance residuals are as defined in Details of \code{\link{llm.fit}} (there is no conceptual distinction between scaled and  unscaled residuals here, since the residual dispersion parameter is not generally a scale factor, but the returned deviance residuals for non-GLMs are analogous to the scaled ones for GLMs as they depend on residual dispersion). \code{"std_dev_res"} and \code{"std_dev_rt"} residuals are defined from them as detailed above for GLM response families, with the additional convention that \eqn{\phi=1} (since the family's own residual dispersion parameter already enters in the definition of deviance residuals for non-GLM families).  Pearson residuals and response residuals are defined as in \code{stats:::residuals.glm}.  The \code{"working"} residuals are defined for each response as \eqn{- [d \log(clik)/d \eta]/[d^2 \log(clik)/d \eta^2]} where clik is the conditional likelihood.
}
\value{ A vector of residuals}
\references{
  Lee, Y., Nelder, J. A. and  Pawitan, Y. (2006). Generalized linear models with random effects: unified analysis via h-likelihood. Chapman & Hall: London.
  
Dunn, K. P., and Smyth, G. K. (1996). Randomized quantile residuals. Journal of Computational and Graphical Statistics 5, 1-10.  
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
