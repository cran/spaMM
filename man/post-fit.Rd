\name{post-fit}
\alias{post-fit}
\alias{RLRsim}
\alias{multcomp}
\alias{glht}
\alias{DHARMa}
\alias{lmerTest}
\title{
Applying post-fit procedures on spaMM results
}
\description{
  Packages implementing post-fit procedures define helper functions which may not handle \pkg{spaMM}'s fit objects, or which have not always handled them, or which can handle them correctly only with some non-default arguments. This documentation topic gives further directions to apply some such post-fit procedures (from packages \pkg{DHARMa}, \pkg{RLRsim}, \pkg{multcomp} and \pkg{lmerTest}) to these fit objects. 
  
  For other procedures not considered here, diagnosing a failure in a debugging session may suggest a simple solution (as it did for \code{multcomp::glht}). 
}
\details{
  For multiple comparison procedures by \verb{multcomp::glht}, one has to explicitly give the argument \code{coef.=fixef.HLfit} (see Examples; \code{fixef.HLfit} is the \pkg{spaMM} method for the generic function \code{fixef});
  
  For \pkg{DHARMa} plots, see Details of \code{\link{plot.HLfit}};
  
  For using \pkg{RLRsim::RLRTSim()}, see \code{\link{get_RLRTSim_args}}.
  
  For using \pkg{lmerTest::contest()} or \pkg{lmerTest::anova()}, see \code{\link{as_LMLT}}.
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