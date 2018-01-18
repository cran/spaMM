\name{HLfit}
\alias{HLfit}
\alias{Beta-distribution-random-effects}
\title{Fit mixed models with given correlation matrix}
\description{
  This function fits GLMMs as well as some hierarchical generalized linear models (HGLM; Lee and Nelder 2001).
  \code{HLfit} fits both fixed effects parameters, and dispersion parameters i.e. the variance of the random effects 
  and the variance of the residual error. The linear predictor is of the 
  standard form \code{offset+ X beta + Z b}, where 
  X is the design matrix of fixed effects and Z is a design matrix of random effects.   
  The function also handles modelling of the residual variance.\cr
  \code{HLfit} uses iterative algorithms that may be slow in some applications, in which case it is worth considering the \code{\link{fitme}} function.
  \code{fitme} and \code{corrHLfit} act by optimizing the likelihood of \code{HLfit} fits for different values of random effect parameters or residual dispersion parameters, with distinct usage of the iterative methods implemented in HLfit: \code{fitme} rests more of generic optimization, while \code{corrHLfit} rests more on iterative methods for estimating variance parameters. 
}
\usage{
HLfit(formula, data, family = gaussian(), rand.family = gaussian(), 
      resid.model = ~1, resid.formula, REMLformula = NULL, 
      verbose = c(trace = FALSE), HLmethod = "HL(1,1)", control.HLfit = list(), 
      control.glm = list(), init.HLfit = list(), ranFix = list(), 
      etaFix = list(), prior.weights = NULL, processed = NULL)
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
  \bold{or} a list, with at most two possible elements if its formula involves only fixed effects:
  \describe{
    \item{formula}{model formula as in formula-only case, without left-hand side}
    \item{family}{Always Gamma, with by default a log link. \code{Gamma(identity)} can be tried but may fail because only the log link ensures that the fitted \eqn{\phi} is positive.}
  }
  and additional possible elements (all named as \code{fitme} arguments) if its formula involves random effects:  see \code{\link{phiHGLM}}.
}
  \item{resid.formula}{
  Obsolete, for back-compatibility; will be deprecated. Same as formula in \code{resid.model}. 
}
  \item{REMLformula}{
  A model \code{formula} that allows the estimation of dispersion parameters, and 
  computation of restricted likelihood (\code{p_bv}) under a model different from the predictor \code{formula}.

  For example, if only random effects are included in \code{REMLformula}, an ML fit is performed and \code{p_bv} equals
  the marginal likelihood (or its approximation), \code{p_v}. This ML fit can be performed more simply by setting 
  \code{HLmethod="ML"} and leaving \code{REMLformula} at its default NULL value.
}
  \item{verbose}{
    A vector of booleans. \code{trace} controls various diagnostic messages (possibly messy, and in bad need of revision) about the iterations. \code{TRACE} is most useful to follow the progress of a long computation, particularly in \code{fitme} or \code{corrHLfit} calls.
  }
  \item{HLmethod}{
  Allowed values are \code{"REML"}, \code{"ML"}, \code{"EQL-"} and \code{"EQL+"} for all models;  
  \code{"PQL"} (=\code{"REPQL"}) and \code{"PQL/L"} for GLMMs only; and further values 
  for those curious to experiment (see Details). \bold{The default is REML} (standard REML for LMMs, 
  an extended definition for other models). REML can be viewed as a fom of conditional inference, and non-standard  conditionings can be called as \dQuote{REML} with a non-standard \code{REMLformula}. See Details for further information.
}
  \item{control.HLfit}{
  A list of parameters controlling the fitting algorithms.

  \code{resid.family} allows one to change the link for modeling of residual variance \eqn{\phi}, which is \code{"log"} by default. The family is always Gamma, so the non-default possible values of \code{resid.family} are \code{Gamma(identity)} or \code{Gamma(inverse)}. Only the default value ensures that the fitted \eqn{\phi} is positive.     

  Controls for the fitting algorithms should be ignored in routine use. They are:

  \code{conv.threshold} and \code{spaMM_tol}: \code{spaMM_tol} is a list of tolerance values, with elements \code{Xtol_rel} and \code{Xtol_abs} that define thresholds for relative and absolute changes in parameter values in iterative algorithms (used in tests of the form \dQuote{d(param)< Xtol_rel * param + Xtol_abs}, so that \code{Xtol_abs} is operative only for small parameter values). \code{conv.threshold} is the older way to control \code{Xtol_rel}. Default values are given by spaMM.getOption("spaMM_tol");
   
   \code{break_conv_logL}, a boolean specifying whether the iterative algorithm should terminate when log-likelihood appears to have converged (roughly, when its relative variation over on iteration is lower than 1e-8). Default is FALSE (convergence is then assessed on the parameter estimates rather than on log-likelihood).
%
  \code{iter.mean.dispFix}, the number of iterations of the iterative algorithm for coefficients of the linear predictor,
       if no dispersion parameters are estimated by the iterative algorithm. Defaults to 200; 
%
  \code{iter.mean.dispVar}, the number of iterations of the iterative algorithm for coefficients of the linear predictor,
       if some dispersion parameter(s) is estimated by the iterative algorithm. Defaults to 50;  
%
  \code{max.iter}, the number of iterations of the iterative algorithm for joint estimation of dispersion parameters and
        of coefficients of the linear predictor. Defaults to 200. This is typically much more than necessary, 
        unless there is little information to separately estimate \eqn{\lambda} and \eqn{\phi} parameters.
}
  \item{control.glm}{
    List of parameters controlling GLM fits, passed to \code{glm.control}; e.g. \code{control.glm=list(maxit=100)}. See \code{\link{glm.control}} for further details.  
  }

  \item{init.HLfit}{
  A list of initial values for the iterative algorithm, with possible elements of the list are 
  \code{fixef} for fixed effect estimates (beta),  
  \code{v_h} for random effects vector \bold{v} in the linear predictor,
  \code{lambda} for the parameter determining the variance of random effects \eqn{u} as drawn from the \code{rand.family} distribution 
  \code{phi} for the residual variance. 
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
   An optional vector of prior weights as in \code{\link{glm}}. This fits the data to a model with residual variance \code{phi/prior.weights}, so that increasing the weights by a constant factor \emph{f} will yield (Intercept) estimates of \code{phi} also increased by \emph{f} (this effect cannot be generally achieved if a non-trivial \code{resid.formula} with log link is used). This is not necessarily the way prior weights are interpreted in widely used packages, but this is consistent with what \code{glm} fits. 
 }
  \item{processed}{
    A list of preprocessed arguments, for programming purposes only (as in \code{corrHLfit}).
 }
}
\details{

 \bold{I. Fitting methods:}
 Many approximations for likelihood have been defined to fit mixed models (e.g. Noh and Lee (2007) for some overview), 
 and this function implements several of them, and some additional ones. In particular, 
 PQL as originally defined by Breslow and Clayton (1993) uses REML to estimate dispersion parameters, but this function allows one to use an ML variant of PQL. Moreover, it allows some non-standard specification of the model formula that determines the conditional distribution used in REML.
 
  EQL stands for the EQL method of Lee and Nelder (2001). The '+' version includes the d v/ d tau correction 
  described p. 997 of that paper, and the '-' version ignores it.
  PQL can be seen as the version of EQL- for GLMMs. It estimates fixed effects by maximizing h-likelihood and dispersion parameters by an approximation of REML, i.e. by maximization of an approximation of restricted likelihood.
  PQL/L is PQL without the leverage corrections that define REML estimation of random-effect parameters. Thus, it estimates dispersion parameters by an approximation of marginal likelihood.

 \code{HLmethod} also accepts values of the form \code{"HL(<...>)"}, \code{"ML(<...>)"} and \code{"RE(<...>)"}, e.g. \code{HLmethod="RE(1,1)"}, which allow a more direct specification of the approximations used. 
 HL and RE are equivalent (both imply an REML correction).   
 The first '1' means that a first order Laplace approximation to the likelihood is used to estimate fixed effects 
 (a '0' would instead mean that the h likelihood is used as the objective function).
  The second  '1' means that a first order Laplace approximation to the likelihood or restricted likelihood 
  is used to estimate dispersion parameters, this approximation including the dv/d tau term specifically discussed by Lee & Nelder 2001, p. 997 (a '0' would instead mean that these terms are ignored).
  
  It is possible to enforce the EQL approximation for estimation of dispersion parameter (i.e., Lee and Nelder's (2001) method) by adding a third index with value 0. \code{"EQL+"} is thus \code{"HL(0,1,0)"}, while \code{"EQL-"} is \code{"HL(0,0,0)"}. \code{"PQL"} is EQL- for GLMMs. \code{"REML"} is \code{"HL(1,1)"}. \code{"ML"} is \code{"ML(1,1)"}. 
  
  Some of these distinctions make sense for \bold{GLMs}, and \code{\link[stats]{glm}} methods use approximations, which make a difference for Gamma GLMs. This means in particular that, (as stated in \code{\link[stats]{logLik}}) the logLik of a Gamma GLM fit by \code{glm} differs from the exact likelihood. Further, the dispersion estimate returned by \code{summary.glm} differs from the one implied by \code{logLik}, because \code{summary.glm} uses Pearson residuals instead of deviance residuals, and no \code{HLmethod} tries to reproduce this behaviour. \code{logLik} gives the approximation returned by an \code{"ML(0,0,0)"} fit. The dispersion estimate returned by an \code{"HL(.,.,0)"} fit matches what can be computed from residual deviance and residual degrees of freedom of a glm fit, but this is not the estimate displayed by \code{summary.glm}. With a log link, the fixed effect estimates are unaffected by these distinctions.     

  \bold{II. Random effects} are constructed in several steps. first, a vector \bold{u} of independent and identically distributed (iid) random effects is drawn from some distribution;
 second, a transformation \bold{v}=f(\bold{u}) is applied to each element (this defines \bold{v} which elements are still iid); third, correlated random effects are obtained as \bold{Lv} 
  where \bold{L} is the \dQuote{square root} of a correlation matrix (this may be meaningful only for Gaussian random effects). Coefficients in a random-coefficient model correspond to \bold{Lv}.  
  Finally, a matrix \bold{Z} (or sometimes \bold{ZA}, see \code{\link{Predictor}}) allows to specify how the correlated random effects
  affect the response values. In particular, \bold{Z} is the identity matrix if there is a single observation (response) for each location, but otherwise
  its elements \eqn{z_{ji}} are 1 for the \eqn{j}th observation in the \eqn{i}th location. 
  The design matrix for \bold{v} is then of the form \bold{ZL}. 

 The specification of the random effects u and v handles the following cases: 
 \bold{Gaussian} with zero mean, unit variance, and identity link; \bold{Beta}-distributed, where \eqn{u ~ B(1/(2\lambda),1/(2\lambda))} with mean=1/2, and var\eqn{=\lambda/[4(1+\lambda)]}; and with logit link \code{v=logit(u)};
 \bold{Gamma}-distributed random effects, where \eqn{u ~ }\code{Gamma(shape=}1+1/\eqn{\lambda},\bold{scale=}1/\eqn{\lambda}): see \code{\link{Gamma}} for allowed links and further details; and \bold{Inverse-Gamma}-distributed random effects, where \eqn{u ~ }inverse-Gamma(\code{shape=}1+1/\eqn{\lambda},\bold{rate=}1/\eqn{\lambda}): see \code{\link{inverse.Gamma}} for allowed links and further details.

  \bold{III. The standard errors} reported may sometimes be misleading. For each set of parameters among \eqn{\beta}, \eqn{\lambda}, and \eqn{\phi} parameters these are computed assuming that the other parameters are known without error. This is why they are labelled \code{Cond. SE} (conditional standard error). This is most uninformative in the unusual case where \eqn{\lambda} and \eqn{\phi} are not separately estimable parameters. Further, the SEs for \eqn{\lambda} and \eqn{\phi} are rough approximations as discussed in particular by Smyth et al. (2001; \eqn{V_1} method).    

}
\value{
An object of class \code{HLfit}, which is a list with many elements, not all of which are documented. 

A few extractor functions are available (see \code{\link{extractors}}), 
and should be used as far as possible as they should be backward-compatible from version 1.4 onwards, while the structure of the return object may still evolve. The following information will be useful for extracting further elements of the object.

Elements include \bold{descriptors of the fit}:

\item{eta}{Fitted values on the linear scale  (including the predicted random effects);}
\item{fv}{Fitted values (\eqn{\mu=}<inverse-link>(\eqn{\eta})) of the response variable (returned by the \code{fitted} function);}
\item{fixef}{The fixed effects coefficients, \eqn{\beta} (returned by the \code{fixef} function);}
\item{ranef}{The random effects \eqn{u} (returned by \code{ranef(*,type="uncorrelated")};}
\item{v_h}{The random effects on the linear scale, \eqn{v};}
\item{phi}{The residual variance \eqn{\phi};}
\item{phi.object}{A possibly more complex object describing \eqn{\phi};}
\item{lambda}{The random effects (\eqn{u}) variance \eqn{\lambda};}
\item{lambda.object}{A possibly more complex object describing \eqn{\lambda};}
\item{corrPars}{Agglomerates information on correlation parameters, either fixed, or estimated by \code{HLfit}, \code{corrHLfit} or \code{fitme};}
\item{APHLs}{A list which elements are various likelihood components, include conditional likelihood, h-likelihood, and the two adjusted profile h-likelihoods: the (approximate) marginal \bold{likelihood} \code{p_v} and the (approximate) \bold{restricted likelihood} \code{p_bv} (the latter two available through the \code{logLik} function). See the extractor function \code{\link{get_any_IC}} for information criteria (\dQuote{AIC}) and effective degrees of freedom;}
\item{}{The covariance matrix of \eqn{\beta} estimates is not included as such, but can be extracted by \code{\link{vcov}};}

\bold{Information about the input} is contained in output elements named as \code{HLfit} or \code{corrHLfit} arguments (\code{data,family,resid.family,ranFix,prior.weights}), with the following notable exceptions or modifications:

\item{predictor}{The linear predictor, including the \code{formula} (possibly reformatted) and several attributes;}
\item{resid.predictor}{Analogous to \code{predictor}, for the residual variance;}
\item{rand.families}{corresponding to the \code{rand.family} input;}

\bold{Further miscellaneous diagnostics and descriptors of model structure:}

\item{X.pv}{The design matrix for fixed effects;}
\item{ZAlist,strucList}{Two lists of matrices, respectively the design matrices \dQuote{\bold{Z}}, and the \dQuote{\bold{L}} matrices, for the different random effects terms. The extractor \code{\link{get_ZALMatrix}} can be used to reconstruct a single \dQuote{\bold{ZL}} matrix for all terms.}
\item{fixef_terms,fixef_levels}{Further information about fixed effect model;}
\item{weights}{(binomial data only) the binomial denominators;}
\item{y}{the response vector; for binomial data, the frequency response.}
\item{models}{Additional information on model structure for \eqn{\eta}, \eqn{\lambda} and \eqn{\phi};}
\item{HL}{A set of indices that characterize the approximations used for likelihood;}
\item{leve_phi,lev_lambda}{Leverages;}
\item{dfs}{degrees of freedom for different components of the model;}
\item{warnings}{A list of warnings for events that may have occurred during the fit.}

Finally, the object includes programming tools: \code{call, spaMM.version, fit_time} and \code{envir}.

}
\references{
Breslow, NE, Clayton, DG. (1993). Approximate Inference in Generalized Linear Mixed Models.
Journal of the American Statistical Association 88, 9-25.

Lee, Y., Nelder, J. A. (2001)  Hierarchical generalised linear models: A
synthesis of generalised linear models, random-effect models and structured
dispersions. Biometrika 88, 987-1006.

Lee, Y., Nelder, J. A. and Pawitan, Y. (2006). Generalized linear models with random effects: unified analysis via
h-likelihood. Chapman & Hall: London.

Noh, M., and Lee, Y. (2007). REML estimation for binary data in GLMMs, J.
Multivariate Anal. 98, 896-915.

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
HLfit(y ~X1+X2+X1*X3+X2*X3+I(X2^2)+(1|batch),family=Gamma(log),
          resid.model = ~ X3+I(X3^2) ,data=wafers)
%- : tested in update.Rd
## Gamma - inverseGamma HGLM with log link
HLfit(y ~X1+X2+X1*X3+X2*X3+I(X2^2)+(1|batch),family=Gamma(log),
          HLmethod="HL(1,1)",rand.family=inverse.Gamma(log),
          resid.model = ~ X3+I(X3^2) ,data=wafers)
}          
\keyword{ model }
