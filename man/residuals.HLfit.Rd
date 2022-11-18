\name{residuals.HLfit}
\alias{residuals}
\alias{residuals.HLfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract model residuals
}
\description{
Extracts several types of residuals from an object of class \code{HLfit}. Note that the default type (\code{"deviance"}) of returned residuals differs from the default  (response residuals) of equivalent functions in base R.  
}
\usage{
\method{residuals}{HLfit}(object, 
  type = c("deviance", "pearson", "response", "working", "std_dev_res"), force=FALSE, ...)
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
The four types \code{"deviance"} (default), \code{"pearson"}, \code{"response"} are \code{"working"} are, for GLM families, the same that are returned by \code{\link[stats]{residuals.glm}}. \code{"working"} residuals may be returned only for fixed-effect models. \code{"deviance"} residuals are the signed square root of those returned by \code{\link{dev_resids}}.

  Following Lee et al. (2006, p.52), the standardized deviance residuals returned for \code{type="std_dev_res"} are defined as the deviance residuals divided by \eqn{\phi\sqrt(1-q)}, where the deviance residuals are defined as for a GLM, \eqn{\phi} is the dispersion parameter of the response family (a vector of values, for heteroscedastic cases), and \eqn{q} is a vector of leverages given by \code{hatvalues(., type="std")} (see \code{\link{hatvalues}} for details about these specific standardizing leverages).
  
  These definitions are special cases of more general ones holding for non-GLM response families. In the latter case, the deviance residuals are as defined in Details of \code{\link{llm.fit}}, and \code{"std_dev_res"} residuals are defined from them as above for GLM response families, with the additional convention that \eqn{\phi=1}.  Pearson residuals and response residuals are defined as in \code{stats:::residuals.glm}.  The \code{"working"} residuals are defined for each response as \eqn{- [d \log(clik)/d \eta]/[d^2 \log(clik)/d \eta^2]} where clik is the conditional likelihood.
}
\value{ A vector of residuals}
\references{
  Lee, Y., Nelder, J. A. and Pawitan, Y. (2006). Generalized linear models with random effects: unified analysis via
  h-likelihood. Chapman & Hall: London.
}
\examples{
data("wafers")
fit <- fitme(y ~X1+(1|batch) ,data=wafers, init=list(phi=NaN))  # : this 'init' 
#                 implies that standardized deviance residuals are saved in the 
#                 fit result, allowing the following comparison: 

r1 <- residuals(fit, type="std_dev_res") # gets stored value
r2 <- residuals(fit, type="std_dev_res", force=TRUE) # forced recomputation
if (diff(range(r1-r2))>1e-14) stop()
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
