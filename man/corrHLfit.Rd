\name{corrHLfit}
\alias{corrHLfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Fits a mixed model, typically a spatial GLMM.
}
\description{
  \code{corrHLfit} performs the joint estimation of correlation parameters, fixed effect and dispersion parameters. 
}
\usage{
corrHLfit(formula, data, init.corrHLfit = list(), 
          init.HLfit = list(), ranFix = list(), lower = list(), 
          upper = list(), 
          objective = NULL, resid.model = ~1, resid.formula,
          control.dist = list(), control.corrHLfit = list(),
          processed = NULL, family = gaussian(), 
          nb_cores = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
Either a linear model \code{\link{formula}} (as handled  by various fitting functions) or a \code{predictor}, i.e. a formula with attributes (see \code{\link{Predictor}} and examples below). See Details in \code{\link{spaMM}} for allowed terms in the formula.
}
  \item{data}{A data frame containing the variables in the response and the model formula.}
  \item{init.corrHLfit}{
An optional list of initial values for correlation and/or dispersion parameters, e.g. 
\code{list(rho=1,nu=1,lambda=1,phi=1)} where \code{rho} and \code{nu} are parameters of the Matérn family (see \code{\link{Matern}}), and 
\code{lambda} and \code{phi} are dispersion parameters (see Details in \code{\link{spaMM}} for the meaning of these parameters). 
 All are optional, but giving values for 
a dispersion parameter changes the ways it is estimated (see Details).
\code{rho}  may be a vector (see \code{\link{make_scaled_dist}}) and, in that case, it is possible that some or all of its elements are \code{NA}, for which \code{corrHLfit} substitute automatically determined values.  
}
  \item{init.HLfit}{See identically named \code{\link{HLfit}} argument.
%- FR->FR more to say (if corr estim within HLfit becomes public ?)...
  }
  \item{ranFix}{A list similar to \code{init.corrHLfit}, but specifying fixed values of the parameters not estimated. See \code{\link{ranFix}} for further information.}
  \item{lower}{
An optional list of values of parameters specified through \code{init.corrHLfit}, used as lower values in calls to \code{optim}. See Details for default values.
}
  \item{upper}{Same as \code{lower}, but upper values.}
  \item{objective}{ Can be ignored: default value is determined from other 
  arguments (\code{HLmethod} and \code{REMLformula}), being \code{"p_bv"} for restricted likelihood or \code{"p_v"} for marginal likelihood. Non-default values can be given but this is an obsolete and the other arguments should be used.}
  \item{resid.model, resid.formula}{See identically named \code{\link{HLfit}} arguments.}
\item{control.dist}{See \code{control.dist} in \code{\link{HLCor}}}
\item{control.corrHLfit}{
  This may be used control the optimizer. See \code{\link{spaMM.options}}} for default values.
  \item{processed}{For programming purposes, not documented.  }
  \item{family}{Either a \code{\link{family}} or a \code{\link{multi}} value. }
  \item{nb_cores}{Number of cores to use for parallel computations.}
  \item{\dots}{
Optional arguments passed to \code{\link{HLCor}}, \code{\link{HLfit}} or  \code{\link{designL.from.Corr}}, for example the \code{distMatrix} argument 
of \code{HLCor}. Arguments that do not fit within these functions are detected and a warning is issued.
}
}
\details{
  Under the Matérn correlation model, corrHLfit typically performs a optimization over the \eqn{\rho} and \eqn{\nu} parameters, 
  with maximum possible values as set by \code{\link{spaMM.options}}.\cr 

  By default \code{corrHLfit} will estimate correlation parameters by maximizing the \code{objective} value returned by \code{HLCor} calls wherein the dispersion parameters are estimated jointly with fixed effects for given correlation parameters. If dispersion parameters are specified in \code{init.corrHLfit}, they will also be estimated by maximizing the \code{objective} value, and \code{HLCor calls} will not estimate them jointly with fixed effects. This means that in general the fixed effect estimates may vary depending on \code{init.corrHLfit} when any form of REML correction is applied. 
  
  Correctly using \code{corrHLfit} for likelihood ratio tests
  of fixed effects may then by tricky. It is safe to perform full ML fits of all parameters (using \code{HLmethod="ML"}) for such tests (see Examples).   The higher level function \code{\link{fixedLRT}} is a safe interface for likelihood ratio tests using some form of REML estimation in \code{corrHLfit}.  
  
  \code{attr(<fitted object>,"optimInfo")$lower} and ...\code{$upper} gives the lower and upper bounds for optimization of correlation parameters. These are the default values if the user did not provide explicit values. For the adjacency model, the default values are the inverse of the maximum and minimum eigenvalues of the \code{adjMatrix}. For the Matérn model, the default values are not so easily summarized: they are intended to cover the range of values for which there is statistical information to distinguish among them.
}
\value{
The return value of an \code{HLCor} call, with additional attributes. The \code{HLCor} call is evaluated at the estimated correlation parameter values. These values are included in the return object as its \code{$corrPars} member. The attributes added by \code{corrHLfit} include the original call of the function (which can be retrived by \code{getCall}(<fitted object>), and information about the optimization call within \code{corrHLfit}. 
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
  See more examples on data set \code{\link{Loaloa}}.
  See \code{\link{fixedLRT}} for likelihood ratio tests.
}
\examples{
# Example with an adjacency matrix (autoregressive model):
# see 'adjacency' documentation page

#### Examples with Matérn correlations
## A likelihood ratio test based on the ML fits of a full and of a null model.
if (spaMM.getOption("example_maxtime")>1.4) {
 data("blackcap")
 fullfit <- corrHLfit(migStatus ~ means+ Matern(1|latitude+longitude),data=blackcap,
                    HLmethod="ML") 
 summary(fullfit)
 nullfit <- corrHLfit(migStatus ~ 1 + Matern(1|latitude+longitude),data=blackcap,
                    HLmethod="ML",init.corrHLfit=list(phi=1e-6)) 
 summary(nullfit)
 ## p-value:
 1-pchisq(2*(logLik(fullfit)-logLik(nullfit)),df=1)
}

## see data set Loaloa for additional examples
}
\keyword{ model }


