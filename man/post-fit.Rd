\name{post-fit}
\alias{post-fit}
\alias{RLRsim}
\alias{multcomp}
\alias{glht}
\alias{DHARMa}
\title{
Applying post-fit procedures on spaMM results
}
\description{
  Packages implementing post-fit procedures define helper functions that may handle a limited range of classes of fit results. This documentation topic gives further directions to apply some post-fit procedures that handle objects of class \code{HLfit}. For other procedures, \code{spaMM} results may or may not be handled correctly by default. In such cases, diagnosing a failure in a debugging session may suggest a simple solution (such as providing the \code{coef.} argument to \verb{multcomp::glht}). 
}
\details{
  For multiple comparison procedures by \verb{multcomp::glht}, one has to explicitly give the argument \code{coef.=fixef.HLfit} (see Examples; \code{fixef.HLfit} is the spaMM method for the generic function \code{fixef});
  
  For \code{DHARMa} plots, see \code{\link{plot.HLfit}};
  
  For using \code{RLRsim::RLRTSim()}, see \code{\link{get_RLRTSim_args}}.
}
\examples{
if (requireNamespace("multcomp", quietly = TRUE)) {
  library(multcomp)
  set.seed(123)
  irisr <- cbind(iris,id=sample(4,replace=TRUE,size=nrow(iris)))
  irisfit <- fitme(Petal.Length~ Species +(1|id), data=irisr, family=Gamma(log))
  summary(glht(irisfit,mcp("Species" = "Tukey"), coef.=fixef.HLfit))
}
}