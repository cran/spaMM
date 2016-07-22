\name{fitme}
\alias{fitme}
\title{
Fitting function for fixed- and mixed-effect models with GLM response. 
}
\description{
This is a common interface for fitting most models that spaMM can fit, therefore substituting to \code{corrHLfit}, \code{HLCor} and \code{HLfit}. Random-slope models are not yet handled by this function. Currently it is fast for some models and slow for others, with few general guidance except that \code{fitme} should generally be faster than the alternative fitting functions for large data sets when the residual variance model is a single constant term (no structured dispersion). 
}
\usage{
fitme(formula, data, family = gaussian(), init = list(), fixed = list(), 
      lower = list(), upper = list(), resid.model = ~1, init.HLfit = list(), 
      control = list(), control.dist = list(), method = "ML", 
      HLmethod = method, processed = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
Either a linear model \code{\link{formula}} (as handled  by various fitting functions) or a \code{predictor}, i.e. a formula with attributes (see \code{\link{Predictor}} and examples below). See Details in \code{\link{spaMM}} for allowed terms in the formula.
}
  \item{data}{
A data frame containing the variables in the response and the model formula.
}
  \item{family}{
  Either a \code{\link{family}} or a \code{\link{multi}} value. 
}
  \item{init}{
An optional list of initial values for correlation and/or dispersion parameters and/or response family parameters, e.g. 
\code{list(rho=1,nu=1,lambda=1,phi=1)} where \code{rho} and \code{nu} are parameters of the Matérn family, and 
\code{lambda} and \code{phi} are dispersion parameters (see Details in \code{\link{spaMM}} for the meaning of these parameters). 
 All are optional, but giving values for a dispersion parameter changes the ways it is estimated (see Details).
\code{rho}  may be a vector (see \code{\link{make_scaled_dist}}) and, in that case, it is possible that some or all of its elements are \code{NA}, for which \code{corrHLfit} substitute automatically determined values.  
}
  \item{fixed}{
A list similar to \code{init}, but specifying fixed values of the parameters not estimated. 
}
  \item{lower}{
An optional list of values of parameters specified through \code{init.corrHLfit}, used as lower values in calls to \code{optim}. See Details for default values.
}
  \item{upper}{Same as \code{lower}, but upper values.}
\item{resid.model}{
See identically named \code{\link{HLfit}} argument.
}
  \item{init.HLfit}{
See identically named \code{\link{HLfit}} argument.
}
  \item{control}{
  A list of control parameters, with one possible element \code{control$nloptr}, itself a list of control parameters to be copied in the \code{opts} argument of \code{\link[nloptr]{nloptr}}. Default controls are\cr 
  \code{list(algorithm="NLOPT_LN_BOBYQA",xtol_rel=1.0e-4,maxeval=-1,print_level=0)}   
}
  \item{control.dist}{
See \code{control.dist} in \code{\link{HLCor}}
}
  \item{method, HLmethod}{
  \code{"ML"} or \code{"REML"}. \code{"ML"} is the default, in contrast to \code{"REML"} for the \code{HLmethod} argument of other fitting functions. Other possible values of \code{HLfit}'s \code{HLmethod} argument are handled and all should give results close to the other fitting methods with the same \code{HLmethod} argument.
}
  \item{processed}{
For programming purposes, not documented.
}
  \item{\dots}{
Optional arguments passed to \code{\link{HLCor}}, \code{\link{HLfit}} or  \code{\link{designL.from.Corr}}, for example the \code{distMatrix} argument of \code{HLCor}. 
% FR->FR : ? :
%Arguments that do not fit within these functions are detected and a warning is issued.
}
}
\details{
\code{fitme} uses \code{\link[nloptr]{nloptr}} for optimization, where \code{corrHLfit} uses \code{optim}.

By default \code{corrHLfit} initiates optimization with the best of a grid of parameter values, while \code{fitme} performs no such grid search for initial value. The grid search may help find a better maximum (as is the case in the toy examples based on the \code{blackcap} data set), but may as well trap the function in an inferior local maximum.  
}
\value{
The return value of an \code{HLCor} or an \code{HLfit} call, with additional attributes. The \code{HLCor} call is evaluated at the estimated correlation parameter values. These values are included in the return object as its \code{$corrPars} member. The attributes added by \code{fitme} include the original call of the function (which can be retrived by \code{getCall}(<fitted object>), and information about the optimization call within \code{fitme}. 
}
\examples{
# see help("COMPoisson")and help("negbin")
}
\keyword{ model }