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
corrHLfit(formula,data,init.corrHLfit=list(),ranFix=list(),lower=list(),upper=list(),
          trace=list(file=NULL,append=T),objective="p_bv",control.corrHLfit=list(),...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
Either a linear model \code{\link{formula}} (as handled  by various fitting functions) or a \code{predictor}, i.e. a formula with attributes (see \code{\link{Predictor}} and examples below).
}
  \item{data}{A data frame containing the variables in the response and the model formula.}
  \item{init.corrHLfit}{
An optional list of initial values for correlation and/or dispersion parameters, e.g. 
\code{list(rho=1,nu=1,lambda=1,phi=1)} where \code{rho} and \code{nu} are parameters of the Matérn family, and 
\code{lambda} and \code{phi} are dispersion parameters (see Details in \code{\link{spaMM}} for the meaning of these parameters). 
 All are optional, but giving values for 
a dispersion parameter changes the ways it is estimated (see Details).
}
  \item{ranFix}{
A list similar to \code{initpar}, but specifying fixed values of the parameters not estimated. 
}
  \item{lower}{
An optional list of values of parameters specified through \code{initpar}, used as lower values in calls to \code{optim}
}
  \item{upper}{
Same as \code{lower}, but upper values.
}
  \item{trace}{
Not for normal use.
If trace is of the form \code{trace=list(file=<filename>,append=F)}, some trace information is written in the file 'filename'. 
This file is written over by each new call of \code{corrHLfit} unless \code{append=T}.
Information is written for each \code{HLCor} call. It includes the \code{APHLs}, and all given fixed random effect 
parameters in the HLCor call. Information about the variables is printed at the end of the file, but may be slightly inaccurate 
(some amount of guesswork is expected from users venturing into \code{trace}). The arguments of all HLCor calls are saved in RData files
(if this is felt inappropriate, then using \code{\link{spaMM.options}(TRACE.UNLINK=TRUE)} will keep only the latest HLCor call).
}
  \item{objective}{
The objective function maximized for estimation of parameters by \code{optim}. 
Either \code{"p_bv"} for restricted likelihood or \code{"p_v"} for marginal likelihood.
}
  \item{control.corrHLfit}{
   Only for development purposes, not documented.
}
  \item{...}{
Optional arguments passed to \code{\link{HLCor}}, \code{\link{HLfit}} or  \code{\link{designL.from.Corr}}, for example the \code{distMatrix} argument 
of \code{HLCor}. Arguments that do not fit within these functions are detected and a warning is issued.
}
}
\details{
  Under the Matérn correlation model, corrHLfit typically performs a optimization over the \eqn{\rho} and \eqn{\nu} parameters, 
  with maximum possible values as set by \code{\link{spaMM.options}}.\cr 

  By default \code{corrHLfit} will estimate correlation parameters by maximizing the \code{objective} value returned by \code{HLCor} calls wherein the dispersion parameters are estimated jointly with fixed effects for given correlation parameters. If dispersion parameters are specified in \code{init.corrHLfit}, they will also be estimated by maximizing the \code{objective} value, and \code{HLCor calls} will not estimate them jointly with fixed effects. This means that in general the fixed effect estimates may vary depending on \code{init.corrHLfit} when any form of REML correction is applied. 
  
  Correctly using \code{corrHLfit} for likelihood ratio tests
  of fixed effects may then by tricky. It is safe to perform full ML fits of all parameters (using \code{objective="p_v",HLmethod="ML"}) for such tests (see Examples).   The higher level function \code{\link{fixedLRT}} is a safe interface for likelihood ratio tests using some form of REML estimation in \code{corrHLfit}.  
}
\value{
A list with currently a single element \code{hlfit}, which is the HLfit object returned by 
\code{HLCor} for the estimated correlation parameter values. These correlation parameter 
values are included in \code{hlfit} as  its \code{$corrPars} member.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
  See more examples on data set \code{\link{Loaloa}}.
  See \code{\link{fixedLRT}} for likelihood ratio tests.
}
\examples{
#### example with an adjacency matrix (autoregressive model)
data(scotlip)
corrHLfit(cases~I(prop.ag/10) +adjacency(1|gridcode)+offset(log(scotlip$expec)),
               data=scotlip,family=poisson(),
               adjMatrix=Nmatrix,lower=list(rho=0),upper=list(rho=0.1745))


#### Examples with Matérn correlations
## A likelihood ratio test based on the ML fits of a full and of a null model.
\dontrun{
data(blackcap)
fullfit <- corrHLfit(migStatus ~ means+ Matern(1|latitude+longitude),data=blackcap,
                    objective="p_v",HLmethod="ML") ## takes ~ 6s
summary(fullfit)
nullfit <- corrHLfit(migStatus ~ 1 + Matern(1|latitude+longitude),data=blackcap,
                    objective="p_v",HLmethod="ML") ## takes ~ 3s
summary(nullfit)
## p-value:
1-pchisq(2*(fullfit$hlfit$APHLs$p_v-nullfit$hlfit$APHLs$p_v),df=1)
}

## see data set Loaloa for additional examples
}
\keyword{ model }


