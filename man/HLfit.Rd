\name{HLfit}
\alias{HLfit}
\alias{REMLformula}
\alias{Beta-distribution-random-effects}
\alias{Beta}
\alias{prior.weights}
\title{Fit mixed models with given correlation matrix}
\description{
  This function fits GLMMs as well as some hierarchical generalized linear models (HGLM; Lee and Nelder 2001).
  \code{HLfit} fits both fixed effects parameters, and dispersion parameters i.e. the variance of the random effects (full covariance for random-coefficient models), and the variance of the residual error. The linear predictor is of the standard form \code{offset+ X beta + Z b}, where X is the design matrix of fixed effects and Z is a design matrix of random effects (typically an incidence matrix with 0s and 1s, but not necessarily). Models are fitted by an iterative algorithm alternating estimation of fixed effects and of dispersion parameters. The residual dispersion may follow a \dQuote{structured-dispersion model} modeling heteroscedasticity. 
  Estimation of the latter parameters is performed by a form of fit of debiased residuals, which allows fitting a structured-dispersion model (Smyth et al. 2001).  However, evaluation of the debiased residuals can be slow in particular for large datasets. For models without structured dispersion, it is then worth using the \code{\link{fitme}} function (or the \code{\link{corrHLfit}} function with non-default arguments). These functions can optimize the likelihood of \code{HLfit} fits for different given values of the dispersion parameters (\dQuote{outer optimization}), thereby avoiding the need to estimate debiased residuals. 
}
\usage{
HLfit(formula, data, family = gaussian(), rand.family = gaussian(), 
      resid.model = ~1, REMLformula = NULL, verbose = c(inner = FALSE), 
      HLmethod = "HL(1,1)", method="REML", control.HLfit = list(), 
      control.glm = list(), init.HLfit = list(), ranFix = list(), 
      etaFix = list(), prior.weights = NULL, weights.form = NULL, processed = NULL)
## see 'rand.family' argument for inverse.Gamma
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
  A \code{\link{formula}}; or a \code{predictor}, i.e. a formula with attributes created by \code{\link{Predictor}}, if design matrices for random effects have to be provided. See Details in \code{\link{spaMM}} for allowed terms in the formula (except spatial ones).
}
  \item{data}{
     A data frame containing the variables named in the model formula.  
}
  \item{family}{
   A \code{family} object describing the distribution of the response variable. See Details in \code{\link{spaMM}} for handled families.
}
  \item{rand.family}{
  A \code{family} object describing the distribution of the random effect, or a \code{list} of 
  family objects for different random effects (see Examples). Possible options are
  \code{gaussian()}, \code{Gamma(log)}, \code{Gamma(identity)} (see Details), \code{Beta(logit)}, \code{inverse.Gamma(-1/mu)}, and \code{inverse.Gamma(log)}.
  For discussion of these alternatives see Lee and Nelder 2001 or Lee et al. 2006, p. 178-.
  Here the family gives the distribution of a random effect \eqn{u} 
  and the link gives \code{v} as function of \eqn{u} (see Details).
  If there are several random effects and only one family is given, this family holds for all random effects.
}
  \item{resid.model}{
  \bold{Either} a formula (without left-hand side) for the dispersion parameter \code{phi} of the residual error. A log link is assumed by default;\cr 
  \bold{or} a list, with at most three possible elements if its formula involves only fixed effects:
  \describe{
    \item{formula}{model formula as in formula-only case, without left-hand side}
    \item{family}{Always Gamma, with by default a log link. \code{Gamma(identity)} can be tried but may fail because only the log link ensures that the fitted \eqn{\phi} is positive.}
    \item{fixed}{can be used to specify the residual dispersion parameter of the residual dispersion model itself. The default value is 1; this argument can be used to set another value, and \code{fixed=list(phi=NA)} will force estimation of this parameter.} 
  }
  and additional possible elements (all named as \code{fitme} arguments) if its formula involves random effects:  see \code{\link{phiHGLM}}.
}
  \item{REMLformula}{
  A model \code{formula} that controls the estimation of dispersion parameters and the computation of restricted likelihood (\code{p_bv}), where the conditioning inherent in REML is defined by a model different from the predictor \code{formula}. A simple example (useless in practice) of its effect is to replicate an ML fit by specifying \code{method="REML"} and an \code{REMLformula} with no fixed effect. The latter implies that no conditioning is performed and that \code{p_bv} equals the marginal likelihood (or its approximation), \code{p_v}. One of the examples in \code{\link{update.HLfit}} shows how \code{REMLformula} can be useful, but otherwise this argument may never be needed for standard REML or ML fits. For non-standard likelihood ratio tests using \code{REMLformula}, see \code{\link{fixedLRT}}.    
}
  \item{verbose}{
    A vector of booleans. The \code{inner} element controls various diagnostic messages (possibly messy) about the iterations. This should be distinguished from the \code{TRACE} element, meaningful in \code{fitme} or \code{corrHLfit} calls, and much more useful. \code{phifit} (which defaults to \code{TRUE}) controls messages about the progress of residual dispersion fits in DHGLMs. 
  }
  \item{method}{Character: the fitting method. 
  allowed values are \code{"REML"}, \code{"ML"}, \code{"EQL-"} and \code{"EQL+"} for all models;  
  \code{"PQL"} (=\code{"REPQL"}) and \code{"PQL/L"} for GLMMs only; and further values 
  for those curious to experiment (see \code{\link{method}}). \bold{The default is REML} (standard REML for LMMs, 
  an extended definition for other models). REML can be viewed as a form of conditional inference, and non-standard conditionings can be called by using a non-standard \code{REMLformula}.
}
  \item{HLmethod}{Same as \code{method}. It is useless to specify \code{HLmethod} when \code{method} is specified. The default value \code{"HL(1,1)"} means the same as \code{method="REML"}, but more accurately relates to definitions of approximations of likelihood in the \eqn{h}-likelihood literature.}
  \item{control.HLfit}{
  A list of parameters controlling the fitting algorithms, which should mostly be ignored in routine use. 
  See \code{\link{control.HLfit}} for possible controls.
}
  \item{control.glm}{
    List of parameters controlling GLM fits, passed to \code{glm.control}; e.g.\cr 
    \code{control.glm=list(maxit=100)}. See \code{\link{glm.control}} for further details.  
  }

  \item{init.HLfit}{
  A list of initial values for the iterative algorithm, with possible elements of the list are 
  \code{fixef} for fixed effect estimates (beta),  
  \code{v_h} for random effects vector \bold{v} in the linear predictor,
  \code{lambda} for the parameter determining the variance of random effects \eqn{u} as drawn from the \code{rand.family} distribution,  
  and \code{phi} for the residual variance. 
  However, this argument can be ignored in routine use. 
}
% F I X M E Break this doc in smaller docs !!!
  \item{ranFix}{
  A list of fixed values of random effect parameters. See \code{\link{ranFix}} for further information.
 }
  \item{etaFix}{
   A list of given values of the coefficients of the linear predictor. See \code{\link{etaFix}} for further information.     
 }
  \item{prior.weights}{
   An optional vector of prior weights as in \code{\link{glm}}. This fits the data to a probability model with residual variance parameter given as \code{phi/prior.weights} instead of the canonical parameter \code{phi} of the response family, and all further outputs are defined to be consistent with this (see section IV in Details). 
 }
  \item{weights.form}{
  Specification of prior weights by a one-sided formula: use \code{weights.form = ~ pw} instead of \code{prior.weights = pw}. The effect will be the same except that such an argument, known to evaluate to an object of class \code{"formula"}, is suitable to enforce safe programming practices (see \code{\link{good-practice}}).  
  }
  \item{processed}{
    A list of preprocessed arguments, for programming purposes only (as in \code{corrHLfit}).
 }
}
\details{

 \bold{I. Approximations of likelihood:} see \code{\link{method}}.

  \bold{II. Possible structure of Random effects:} see \code{\link{random-effects}}, but note that \code{HLfit} does not fit models with autocorrelated random effects.

  \bold{III. The standard errors} reported may sometimes be misleading. For each set of parameters among \eqn{\beta}, \eqn{\lambda}, and \eqn{\phi} parameters these are computed assuming that the other parameters are known without error. This is why they are labelled \code{Cond. SE} (conditional standard error). This is most uninformative in the unusual case where \eqn{\lambda} and \eqn{\phi} are not separately estimable parameters. Further, the SEs for \eqn{\lambda} and \eqn{\phi} are rough approximations as discussed in particular by Smyth et al. (2001; \eqn{V_1} method).    
  
  \bold{IV. prior weights}. This controls the likelihood analysis of heteroscedastic models. In particular, changing the weights by a constant factor \emph{f} should, and will, yield a fit with unchanged likelihood and (Intercept) estimates of \code{phi} also increased by \emph{f} (except if a non-trivial \code{resid.formula} with log link is used). This is consistent with what \code{glm} does, but other packages may not follow this logic (whatever their documentation may say: check by yourself by changing the weights by a constant factor).

}
\value{
An object of class \code{HLfit}, which is a list with many elements, not all of which are documented. 

A few extractor functions are available (see \code{\link{extractors}}), 
and should be used as far as possible as they should be backward-compatible from version 1.4 onwards, while the structure of the return object may still evolve. The following information will be useful for extracting further elements of the object.

Elements include \bold{descriptors of the fit}:

\item{eta}{Fitted values on the linear scale  (including the predicted random effects). \code{predict(.,type="link")} can be used as a formal extractor;}
\item{fv}{Fitted values (\eqn{\mu=}<inverse-link>(\eqn{\eta})) of the response variable. \code{fitted(.)} or \code{predict(.)} can be used as formal extractors;}
\item{fixef}{The fixed effects coefficients, \eqn{\beta} (returned by the \code{fixef} function);}
\item{v_h}{The random effects on the linear scale, \eqn{v}, with atttribute the random effects \eqn{u} (returned by \code{ranef(*,type="uncorrelated")};}
\item{phi}{The residual variance \eqn{\phi};}
\item{phi.object}{A possibly more complex object describing \eqn{\phi};}
\item{lambda}{The random-effect (\eqn{u}) variance(s) \eqn{\lambda} in compact form;}
\item{lambda.object}{A possibly more complex object describing \eqn{\lambda};}
\item{ranef_info}{environment where information about the structure of random effects is stored;}
\item{corrPars}{Agglomerates information on correlation parameters, either fixed, or estimated by \code{HLfit}, \code{corrHLfit} or \code{fitme};}
\item{APHLs}{A list which elements are various likelihood components, include conditional likelihood, h-likelihood, and the Laplace approximations: the (approximate) marginal \bold{likelihood} \code{p_v} and the (approximate) \bold{restricted likelihood} \code{p_bv} (the latter two available through the \code{logLik} function). See the extractor function \code{\link{get_any_IC}} for information criteria (\dQuote{AIC}) and effective degrees of freedom;}

The covariance matrix of \eqn{\beta} estimates is not included as such, but can be extracted by \code{\link{vcov}}.

\bold{Information about the input} is contained in output elements named as \code{HLfit} or \code{corrHLfit} arguments (\code{data,family,resid.family,ranFix,prior.weights}), with the following notable exceptions or modifications:

\item{predictor}{The \code{formula}, possibly reformatted;}
\item{resid.predictor}{Analogous to \code{predictor}, for the residual variance;}
\item{rand.families}{corresponding to the \code{rand.family} input;}

\bold{Further miscellaneous diagnostics and descriptors of model structure:}

\item{X.pv}{The design matrix for fixed effects;}
\item{ZAlist,strucList}{Two lists of matrices, respectively the design matrices \dQuote{\bold{Z}}, and the \dQuote{\bold{L}} matrices, for the different random-effect terms. The extractor \code{\link{get_ZALMatrix}} can be used to reconstruct a single \dQuote{\bold{ZL}} matrix for all terms.}
\item{BinomialDen}{(binomial data only) the binomial denominators;}
\item{y}{the response vector; for binomial data, the frequency response.}
\item{models}{Additional information on model structure for \eqn{\eta}, \eqn{\lambda} and \eqn{\phi};}
\item{HL}{A set of indices that characterize the approximations used for likelihood;}
\item{leve_phi,lev_lambda}{Leverages;}
\item{dfs}{list (possibly structured): degrees of freedom for different components of the model;}
\item{how}{A list containing the information properly extracted by the \code{\link{how}} function.}
\item{warnings}{A list of warnings for events that may have occurred during the fit.}

Finally, the object includes programming tools: \code{call, spaMM.version, fit_time} and \code{envir}.

}
\references{
Lee, Y., Nelder, J. A. (2001)  Hierarchical generalised linear models: A
synthesis of generalised linear models, random-effect models and structured
dispersions. Biometrika 88, 987-1006.

Lee, Y., Nelder, J. A. and Pawitan, Y. (2006). Generalized linear models with random effects: unified analysis via
h-likelihood. Chapman & Hall: London.

Smyth GK, Huele AF, Verbyla AP (2001). Exact and approximate REML for heteroscedastic regression. Statistical Modelling 1, 161-175. 
}
\seealso{
\code{\link{HLCor}} for estimation with given spatial correlation parameters;
\code{\link{corrHLfit}} for joint estimation with spatial correlation parameters;
\code{\link{fitme}} as an alternative to all these functions.
}

\examples{
data("wafers")
## Gamma GLMM with log link
% example also in main page...
HLfit(y ~ X1+X2+X1*X3+X2*X3+I(X2^2)+(1|batch), family=Gamma(log),
          resid.model = ~ X3+I(X3^2) ,data=wafers)
%- : tested in update.Rd
## Gamma - inverseGamma HGLM with log link
HLfit(y ~ X1+X2+X1*X3+X2*X3+I(X2^2)+(1|batch), family=Gamma(log),
          rand.family=inverse.Gamma(log),
          resid.model = ~ X3+I(X3^2) , data=wafers)
}          
\keyword{ model }
