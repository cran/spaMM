\name{fitme}
\alias{fitme}
\title{
Fitting function for fixed- and mixed-effect models with GLM response. 
}
\description{
This is a common interface for fitting most models that spaMM can fit, from linear models to mixed models with non-gaussian random effects, therefore substituting to \code{corrHLfit}, \code{HLCor} and \code{HLfit}. By default, it uses ML rather than REML (differing in this respect from the other fitting functions). It may use \dQuote{outer optimization}, i.e., generic optimization methods for estimating all dispersion parameters, rather than the often slow default method in \code{HLfit}. Outer optimization should generally be faster than the alternative algorithms for large data sets when the residual variance model is a single constant term (no structured dispersion). Which method \code{fitme} actually uses depends on its input arguments, and decision criteria are subject to change in the future. For mixed models, \code{fitme} by default tries to select the fastest method when both can be applied, but may not always succeed in this respect.  
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
\code{list(rho=1,nu=1,lambda=1,phi=1)} where \code{rho} and \code{nu} are parameters of the Matérn family (see \code{\link{Matern}}), and 
\code{lambda} and \code{phi} are dispersion parameters (see Details in \code{\link{spaMM}} for the meaning of these parameters). 
 All are optional, but giving values for a dispersion parameter changes the ways it is estimated (see Details and Examples).
\code{rho}  may be a vector (see \code{\link{make_scaled_dist}}) and, in that case, it is possible that some or all of its elements are \code{NA}, for which \code{corrHLfit} substitute automatically determined values.  
}
  \item{fixed}{
A list similar to \code{init}, but specifying fixed values of the parameters not estimated. See \code{\link{fixed}} for further information.
}
  \item{lower}{
An optional list of values of parameters specified through \code{init.corrHLfit}, used as lower values in calls to \code{optim}. See Details for default values.
}
  \item{upper}{Same as \code{lower}, but upper values.}
\item{resid.model}{ See identically named \code{\link{HLfit}} argument. }
  \item{init.HLfit}{ See identically named \code{\link{HLfit}} argument. }
  \item{control}{   A list of control parameters, with two possible elements: 
  \itemize{
    \item \code{$nloptr}, itself a list of control parameters to be copied in the \code{opts} argument of \code{\link[nloptr]{nloptr}}. Default controls are given by \code{spaMM.getOption('nloptr')}
    \item \code{$refit}, a boolean, or a list of booleans with possible elements \code{$phi} and \code{$lambda}. If either element is set to TRUE, then the corresponding parameters are refitted by the internal \code{HLfit} methods (see Details). If \code{$refit} is TRUE, all \code{phi} and \code{lambda} parameters are refit. By default none of them is refitted, but this default may change in the future. 
  }   
}
  \item{control.dist}{ See \code{control.dist} in \code{\link{HLCor}} }
  \item{method, HLmethod}{
  \code{"ML"} or \code{"REML"}. \code{"ML"} is the default, in contrast to \code{"REML"} for the \code{HLmethod} argument of other fitting functions. Other possible values of \code{HLfit}'s \code{HLmethod} argument are handled and should give results close to the other fitting methods with the same \code{HLmethod} argument, the most notable exception being \code{PQL/L} for binary response data, as \code{fitme} does not take into account the specific way leverages are computed in \code{PQL/L}.
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
For \code{phi} and \code{lambda}, \code{fitme} may not use the internal fitting methods of \code{HLfit}. The latter methods are well suited for structured dispersion models, but require the computation of the so-called leverages, which can be slow. However, they also provide some more information such as the \dQuote{cond. SE} (about which see warning in Details of \code{\link{HLfit}}). To force the evaluation of such information from a \code{fitme} call, use the \code{control$refit} argument (see Example).    
}
\value{
The return value of an \code{HLCor} or an \code{HLfit} call, with additional attributes. The \code{HLCor} call is evaluated at the estimated correlation parameter values. These values are included in the return object as its \code{$corrPars} member. The attributes added by \code{fitme} include the original call of the function (which can be retrived by \code{getCall}(<fitted object>), and information about the optimization call within \code{fitme}. 
}
\examples{
## Contrasting different optimization methods:
# We simulate Gamma deviates with mean mu=3 and variance=2, 
#  ie. phi= var/mu^2= 2/9 in the (mu, phi) parametrization of a Gamma 
#  GLM; and shape=9/2, scale=2/3 in the parametrisation of rgamma().
#  Note that phi is not equivalent to scale: 
#  shape = 1/phi and scale = mu*phi.
set.seed(123)
gr <- data.frame(y=rgamma(100,shape=9/2,scale=2/3))
# Here fitme uses HLfit methods which provide cond. SE for phi by default:
fitme(y~1,data=gr,family=Gamma(log))
# To force outer optimization of phi, use the init argument:
fitme(y~1,data=gr,family=Gamma(log),init=list(phi=1))
# To obtain cond. SE for phi after outer optimization, use the 'refit' control:
fitme(y~1,data=gr,family=Gamma(log),,init=list(phi=1),
      control=list(refit=list(phi=TRUE))) ## or ...refit=TRUE...

## see help("COMPoisson"), help("negbin"), help("Loaloa"), etc., for further examples.
}
\keyword{ model }