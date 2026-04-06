\name{fitme}
\alias{fitme}
\title{
Fitting function for fixed- and mixed-effect models with GLM response. 
}
\description{
This is a common interface for fitting most models that \pkg{spaMM} can fit, from linear models to mixed models with non-gaussian random effects, except multivariate-response models (see \code{\link{fitmv}} for the latter). By default, it uses ML rather than REML (differing in this respect from the other fitting functions \code{corrHLfit}, \code{HLCor} and \code{HLfit}). It tends to use \dQuote{outer optimization}, i.e., generic optimization methods for estimating all dispersion parameters, rather than the iterative methods implemented in \code{HLfit} (see Details). 
}
\usage{
fitme(formula, data, family = gaussian(), init = list(), fixed = list(), 
      lower = list(), upper = list(), resid.model = ~1, init.HLfit = list(), 
      control = list(), control.dist = list(), method = "ML", 
      HLmethod = method, processed = NULL, nb_cores = NULL, objective = NULL, 
      weights.form = NULL, ...)
# '...' may notably include arguments 
#     'rand.family', 'control.HLfit', and 'verbose' 
# (see HLfit() documentation for them), 
#     'covStruct', 'corrMatrix', 'adjMatrix' and 'distMatrix' 
# (see HLCor() documentation for them).
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
An optional list of initial values for correlation and/or dispersion parameters for distribution family. The general syntax is the same as for the \code{\link{fixed}} argument, and a simplified syntax may be used when it is not ambiguous. E.g.,
\code{list(rho=1,nu=1,lambda=1,phi=1)} can be used to fix the \code{rho} and \code{nu} parameters of the Matérn family (see \code{\link{Matern}}), as well as the dispersion parameters denoted
\code{lambda} and \code{phi} in \code{\link{spaMM}} (see Examples below, and in the additional \code{\link{inits}} documentation). 
 All initial values are optional, but giving values for a dispersion parameter changes the ways it is estimated (see Details and Examples).
\code{rho}  may be a vector (see \code{\link{make_scaled_dist}}) and, in that case, it is possible that some or all of its elements are \code{NA}, for which \code{fitme} substitutes automatically determined initial values.  
}
  \item{fixed}{
A list similar to \code{init}, but specifying fixed values of the parameters not estimated. See \code{\link{fixed}} for further information; and keep in mind that fixed fixed-effect coefficients can be passed as the \code{etaFix} argument as part of the \sQuote{\ldots}.
}
  \item{lower, upper}{
optional (sub)lists of values of parameters, in the same format as  \code{init}, used to control lower or upper values in calls to \code{nloptr} or similar generic numerical optimization functions. See \code{\link{optimBounds}} for further information.
}
\item{resid.model}{ See identically named \code{\link{HLfit}} argument. }
  \item{init.HLfit}{ See identically named \code{\link{HLfit}} argument. }
  \item{control.dist}{ See \code{control.dist} in \code{\link{HLCor}} }
  \item{method, HLmethod}{Character: the fitting method to be used, such as \code{"ML"}, \code{"REML"} or \code{"PQL/L"}. \code{"ML"} is the default, in contrast to \code{"REML"} for \code{HLfit}, \code{HLCor} and \code{corrHLfit}. Other possible values of \code{HLfit}'s \code{\link{method}} argument are handled. \code{method=c(<"ML" or "REML">,"exp")} can be distinctly useful for slow fits of models with \code{Gamma(log)} family (see \code{\link{method}}). 
}
  \item{weights.form}{
  Specification of prior weights by a one-sided formula: use \code{weights.form = ~ pw} instead of \code{prior.weights = pw}. The effect will be the same except that such an argument, known to evaluate to an object of class \code{"formula"}, is suitable to enforce safe programming practices (see \code{\link{good-practice}}).  
 }
   \item{control}{   A list of (rarely needed) control parameters, with possible elements: 
  \itemize{
    \item \code{refit}, a boolean, or a \code{list} of booleans with possible elements \code{phi}, \code{lambda} and \code{ranCoefs}. If either element is set to TRUE, then the corresponding parameters are refitted by the internal \code{HLfit} methods (see Details), unless these methods were already selected for such parameters in the main fitting step. If \code{refit} is a single boolean, it affects all parameters. By default no parameter is refitted. 
    \item \code{optimizer}, the numerical optimizer, specified as a string and whose default is controlled by the global \pkg{spaMM} option \code{"optimizer"}. Possible values are \code{"nloptr"}, \code{"bobyqa"}, \code{"L-BFGS-B"} and \code{".safe_opt"}, whose meanings are detailed in the documentation for the \code{optimizer} argument of \code{\link{spaMM.options}}. Better left unchanged unless suspect fits are obtained. 
    \item \code{nloptr}, itself a list of control parameters to be copied in the \code{opts} argument of \code{\link[nloptr]{nloptr}}. Default value is given by \code{spaMM.getOption('nloptr')} and possibly other global \pkg{spaMM} options. Better left unchanged unless you are ready to inspect source code.
    \item \code{bobyqa, optim}, lists of controls similar to \code{nloptr} but for methods \code{"bobyqa"} and \code{"L-BFGS-B"}, respectively.
  }   
}
  \item{nb_cores}{ For development purpose, not documented.} % that's for probitgem: Number of cores to use for parallel computations
  \item{processed}{ For programming purpose, not documented.}
  \item{objective}{ For development purpose, not documented.}
  \item{\dots}{
Optional arguments, notably those listed in the comments of the Usage section above, passed to (or operating as if passed to) either \code{\link{HLfit}} or  \code{\link{HLCor}} or  \code{\link{mat_sqrt}}. 
}
}
\details{
 For approximations of likelihood, see \code{\link{method}}. For the possible structures of random effects, see \code{\link{random-effects}},

For \code{phi}, \code{lambda}, and \code{ranCoefs}, \code{fitme} may or may not use the internal iterative fitting methods of \code{HLfit}. The latter methods are well suited for heteroscedastic (or \dQuote{structured-dispersion}) models, but require computations which can be slow for large datasets. For mixed models, \code{fitme} by default tries to select the fastest method when both can be applied. Therefore, \code{fitme} tends to use generic \dQuote{outer} optimization methods, rather than \dQuote{inner} iterative algorithms, by default for large datasets, unless there is a non-trivial \code{resid.model}. The precise criteria for selection of either method are liable to future changes.  The alternative function \code{corrHLfit} tends to use iterative methods as much as possible.

Further, the internal fitting methods of \code{HLfit} also provide more information such as the \dQuote{cond. SE} (about which see warning in Details of \code{\link{HLfit}}). To force the evaluation of such information after an outer-optimization by a \code{fitme} call, use the \code{control$refit} argument (see Example). Alternatively (and possibly of limited use), one can force inner-optimization of \code{lambda} for a given random effect, or of \code{phi}, by setting it to \code{NaN} in \code{init} (see Example using \sQuote{blackcap} data). The same syntax may be tried for \code{phi}.

The results of REML fits of non-gaussian mixed models by outer-optimization may (generally slightly) differ from those by the iterative method. Explicit values of the \code{init} argument may be used to control the selection of either method in \code{fitme} or \code{corrHLfit}, and may be used to ensure better consistency over successive versions of \code{spaMM}: use \code{NaN} as noted above to enforce inner-optimization, and \code{NA} or a numeric value to enforce outer-optimization.  
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