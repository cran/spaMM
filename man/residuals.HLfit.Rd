\name{residuals.HLfit}
\alias{residuals}
\alias{residuals.HLfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract model residuals
}
\description{
Extracts several types of residuals from an object of class \code{HLfit}.
}
\usage{
\method{residuals}{HLfit}(object, 
  type = c("deviance", "pearson", "response", "std_dev_res"), force=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of class \code{HLfit}, as returned by the fitting functions in \code{spaMM}.}
  \item{type}{
  The type of residuals which should be returned. The alternatives are: \code{"deviance"} (default), \code{"pearson"}, \code{"response"} and \code{"std_dev_res"}. See \code{\link{residuals.glm}} for more information about the first three, and the Details for the last one.
  }
  \item{force}{Boolean: to force recomputation of the \code{"std_dev_res"} residuals even if they are available in the object, for checking purposes.}
  \item{\dots}{For consistency with the generic.}
}
\details{
  Following Lee et al. (2006, p.52), the standardized deviance residuals returned for \code{type="std_dev_res"} are defined as the deviance residuals divided by \eqn{\phi\sqrt(1-q)}, where the deviance residuals are defined as for a GLM, \eqn{\phi} is the dispersion parameter of the response family (a vector of values, for heteroscedastic cases), and \eqn{q} is a vector of leverages given by \code{hatvalues(., type="std")} (see \code{\link{hatvalues}} for details about these specific standardizing leverages).
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
