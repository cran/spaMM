\name{fitme}
\alias{fitme}
\title{
Fitting function for fixed- and mixed-effect models with GLM response. 
}
\description{
This is a common interface for fitting most models that spaMM can fit, from linear models to mixed models with non-gaussian random effects, therefore substituting to \code{corrHLfit}, \code{HLCor} and \code{HLfit}. By default, it uses ML rather than REML (differing in this respect from the other fitting functions). It may use \dQuote{outer optimization}, i.e., generic optimization methods for estimating all dispersion parameters, rather than the iterative methods implemented in \code{HLfit}. The results of REML fits of non-gaussian mixed models by these different methods may (generally slightly) differ. Outer optimization should generally be faster than the alternative algorithms for large data sets when the residual variance model is a single constant term (no structured dispersion). For mixed models, \code{fitme} by default tries to select the fastest method when both can be applied, but precise decision criteria are subject to change in the future. \code{corrHLfit} (with non-default arguments to control the optimization method most suitable to a particular problem) may be used to ensure better consistency over successive versions of \code{spaMM}.
}
\usage{
fitme(formula, data, family = gaussian(), init = list(), fixed = list(), 
      lower = list(), upper = list(), resid.model = ~1, init.HLfit = list(), 
      control = list(), control.dist = list(), method = "ML", 
      HLmethod = method, processed = NULL, nb_cores = NULL, objective = NULL, 
      ...)
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
  Either a response \code{\link{family}} or a \code{\link{multi}} value. 
}
  \item{init}{
An optional list of initial values for correlation and/or dispersion parameters and/or response family parameters, e.g. 
\code{list(rho=1,nu=1,lambda=1,phi=1)} where \code{rho} and \code{nu} are parameters of the Matérn family (see \code{\link{Matern}}), and 
\code{lambda} and \code{phi} are dispersion parameters (see Details in \code{\link{spaMM}} for the meaning of these parameters). 
 All are optional, but giving values for a dispersion parameter changes the ways it is estimated (see Details and Examples).
\code{rho}  may be a vector (see \code{\link{make_scaled_dist}}) and, in that case, it is possible that some or all of its elements are \code{NA}, for which \code{fitme} substitutes automatically determined values.  
}
  \item{fixed}{
A list similar to \code{init}, but specifying fixed values of the parameters not estimated. See \code{\link{fixed}} for further information; and keep in mind that fixed fixed-effect coefficients can be passed as the \code{etaFix} argument as part of the \sQuote{\ldots}.
}
  \item{lower}{
An optional (sub)list of values of the parameters specified through \code{init}, in the same format as  \code{init}, used as lower values in calls to \code{optim}. See Details for default values.
}
  \item{upper}{Same as \code{lower}, but for upper values.}
\item{resid.model}{ See identically named \code{\link{HLfit}} argument. }
  \item{init.HLfit}{ See identically named \code{\link{HLfit}} argument. }
  \item{control}{   A list of control parameters, with two possible elements: 
  \itemize{
    \item \code{nloptr}, itself a list of control parameters to be copied in the \code{opts} argument of \code{\link[nloptr]{nloptr}}. Default controls are given by \code{spaMM.getOption('nloptr')}
    \item \code{refit}, a boolean, or a list of booleans with possible elements \code{$phi}, \code{$lambda} and \code{$ranCoefs}. If either element is set to TRUE, then the corresponding parameters are refitted by the internal \code{HLfit} methods (see Details). If \code{$refit} is a single boolean, it affects of parameters. By default only \code{lambda} is refitted, but this default may change in the future. 
  }   
}
  \item{control.dist}{ See \code{control.dist} in \code{\link{HLCor}} }
  \item{method, HLmethod}{Character: the fitting method to be used, such as \code{"ML"}, \code{"REML"} or \code{"PQL/L"}. \code{"ML"} is the default, in contrast to \code{"REML"} for \code{HLfit}, \code{HLCor} and \code{corrHLfit}. Other possible values of \code{HLfit}'s \code{method} argument are handled.
}
  \item{nb_cores}{\bold{Not yet operative}, only for development purposes. Number of cores to use for parallel computations.}
  \item{processed}{ For programming purposes, not documented.}
  \item{objective}{ For development purpose, not documented.}
  \item{\dots}{
Optional arguments passed to (or operating as if passed to) \code{\link{HLCor}}, \code{\link{HLfit}} or  \code{\link{mat_sqrt}}, for example \code{rand.family}, \code{control.HLfit} or the \code{distMatrix} argument of \code{HLCor} (so that estimation of Matern or Cauchy parameters can be combined with use of an ad hoc distance matrix). 
% FIXME : ? :
%Arguments that do not fit within these functions are detected and a warning is issued.
}
}
\details{
 For approximations of likelihood, see \code{\link{method}}. For the possible structures of random effects, see \code{\link{random-effects}},

For \code{phi}, \code{lambda}, and \code{ranCoefs}, \code{fitme} may or may not use the internal fitting methods of \code{HLfit}. The latter methods are well suited for structured dispersion models, but require computations which can be slow for large datasets. Therefore, \code{fitme} tends to outer-optimize by default for large datasets, unless there is a non-trivial \code{resid.model}. The precise criteria for selection of default method by \code{fitme} are liable to future changes. 

Further, the internal fitting methods of \code{HLfit} also provide some more information such as the \dQuote{cond. SE} (about which see warning in Details of \code{\link{HLfit}}). To force the evaluation of such information after an outer-optimization by a \code{fitme} call, use the \code{control$refit} argument (see Example). Alternatively (and possibly of limited use), one can force inner-optimization of \code{lambda} for a given random effect, or of \code{phi}, by setting it to \code{NaN} in \code{init} (see Example using \sQuote{blackcap} data). The same syntax may be tried for \code{phi}.
}
\value{
The return value of an \code{HLCor} or an \code{HLfit} call, with additional attributes. The \code{HLCor} call is evaluated at the estimated correlation parameter values. These values are included in the return object as its \code{$corrPars} member. The attributes added by \code{fitme} include the original call of the function (which can be retrived by \code{getCall}(<fitted object>), and information about the optimization call within \code{fitme}. 
}
\examples{
## Examples with Matern correlations
## A likelihood ratio test based on the ML fits of a full and of a null model.
 data("blackcap")
 (fullfit <- fitme(migStatus ~ means+ Matern(1|longitude+latitude),data=blackcap) )
 (nullfit <- fitme(migStatus ~ 1 + Matern(1|longitude+latitude),data=blackcap)) 
 ## p-value:
 1-pchisq(2*(logLik(fullfit)-logLik(nullfit)),df=1)

## See ?spaMM for examples of conditional autoregressive model and of non-spatial models. 

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

## Outer-optimization is not necessarily the best way to find a global maximum, 
#  particularly when there is little statistical information in the data:  
if (spaMM.getOption("example_maxtime")>1.6) {
  data("blackcap")
  fitme(migStatus ~ means+ Matern(1|longitude+latitude),data=blackcap) # poor
  #  Compare with the following two ways of avoiding outer-optimization of lambda:
  corrHLfit(migStatus ~ means+ Matern(1|longitude+latitude),data=blackcap,
            method="ML")
  fitme(migStatus ~ means+ Matern(1|longitude+latitude),data=blackcap, 
        init=list(lambda=NaN))
}

## see help("COMPoisson"), help("negbin"), help("Loaloa"), etc., for further examples.
}
\keyword{ model }