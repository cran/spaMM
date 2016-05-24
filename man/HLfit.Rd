\name{HLfit}
\alias{HLfit}
\alias{inverse.Gamma}
\title{Fit mixed models with given correlation matrix}
\description{
  This fonction fits GLMMs as well as some hierarchical generalized linear models (HGLM; Lee and Nelder 2001).
  HLfit fits both fixed effects parameters, and dispersion parameters i.e. the variance of the random effects 
  and the variance of the residual error. The linear predictor is of the 
  standard form \code{offset+ X beta + Z v}, where 
  X is the design matrix of fixed effects and Z is a design matrix of random effects.   
  The function also handles a linear predictor (with only fixed effects) for the residual variance. 
}
\usage{
HLfit(formula, data, family = gaussian(), rand.family = gaussian(), 
      resid.model = ~1, resid.formula, REMLformula = NULL, 
      verbose = c(warn = TRUE, trace = FALSE, summary = FALSE), 
      HLmethod = "HL(1,1)", control.HLfit = list(), 
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
   A \code{family} object describing the distribution of the response variable. 
   Possible values include the \code{gaussian}, \code{poisson}, \code{binomial}, \code{Gamma}, \code{negbin} and 
   \code{COMPoisson} families. Possible combinations of family and link are those allowed by each of these families
   (see \code{\link{family}} for the first four, and specific documentation pages for the last two).
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
  \bold{Either} a formula (without left-hand side) for the variance \code{phi} of the residual error, \bold{or} a list 
  with such a formula as element \code{formula}, and the optional element \code{family} for the model family of the residual model. 
  Currently the formula can only contain fixed effects, including an offset. A log link is assumed by default, but an identity link can be requested by setting \code{resid.model$family} to \code{Gamma(identity)}. 
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
    A vector of booleans. \code{trace} controls various diagnostic (possibly messy) messages about the iterations.
    \code{summary} controls whether a summary of the fit is called by \code{HLfit}.
    \code{warn} is for programming purposes and best ignored.
  }
  \item{HLmethod}{
  Allowed values are \code{"REML"}, \code{"ML"}, \code{"EQL-"} and \code{"EQL+"} for all models;  
  \code{"PQL"} (=\code{"REPQL"}) and \code{"PQL/L"} for GLMMs only;
  and (only for those curious to experiment) expressions of the 
  form \code{"HL(<...>)"}, \code{"ML(<...>)"} and \code{"RE(<...>)"}. HL and RE are equivalent.   
  The default behaviour is RE(1,1), which \bold{by default performs REML} (standard REML for LMMs, 
  an extended definition for other models). But it can also perform non-standard forms of REML 
  (indeed including ML), depending on the \code{REMLformula} given.

  EQL stands for the EQL method of Lee and Nelder (2001). The '+' version includes the d v/ d tau correction 
  described p. 997 of that paper, and the '-' version ignores it.
  PQL can be seen as the version of EQL- for GLMMs.
  PQL/L is PQL without the leverage corrections that define REML estimation of random-effect parameters.
  
  For \bold{GLMs}, the default is still REML. For binomial and Poisson GLMs, this cannot be distinguished from an exact ML fit. 
  Note that the \code{\link{glm}} function performs an EQL analysis, which will differ from the REML (and ML) one for Gamma GLMs.  
  
  See Details for the more general syntax.
}
  \item{control.HLfit}{
  A list of parameters controlling the fitting algorithms.

  \code{resid.family} allows one to change the link for modeling of residual variance \eqn{\phi}, which is \code{"log"} by default. The family is always Gamma, so the non-default possible values of \code{resid.family} are \code{Gamma(identity)} or \code{Gamma(inverse)}. Only the default value ensures that the fitted \eqn{\phi} is positive.     

  Controls for the fitting algorithms should be ignored in routine use. They are
%
  \code{conv.threshold}, 
   a convergence threshold for the iterative algorithm, which controls whether 
   linear predictor terms (fixed effects and inferred random effects), and dispersion parameter estimates have converged. Defaults to 1e-05;
   
   \code{break_conv_logL}, a boolean specifying whether the iterative algorithm should terminate when log-likelihood appears to have convergence (roughly, when its relative variation over on iteration is lower than 1e-8). Default is FALSE.
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
  \item{ranFix}{
  A list of fixed values of random effect parameters, with possible elements \code{lambda}, and also \code{phi} for gaussian and Gamma HGLMs. 
  Inhibits the estimation of these parameters. 
 }
  \item{etaFix}{
   A list of fixed values of the coefficients of the linear predictor, with currently documented element \code{beta}. \code{etaFix$beta} should be a vector with names matching (a subset of) coefficient names of a fit without fixed values. It provides a convenient interface for fixing (some of) the fixed-effects coefficients (\eqn{\beta}). In contrast to an offset specification, it affects the REML correction for estimation of dispersion parameters, which depends only on which \eqn{\beta} coefficients are estimated. However, for non-standard use, REML can still be performed as if all \eqn{\beta} coefficients were estimated, by adding attribute \code{keepInREML=TRUE} to \code{etaFix$beta} (see Examples). These different behaviours will be overridden whenever a non-null \code{REMLformula} is provided.     
 }
  \item{prior.weights}{
   An optional vector of prior weights as in \code{\link{glm}}. This fits the data to a model with residual variance \code{phi/prior.weights}, so that increasing the    weights by a constant factor \emph{f} will yield (Intercept) estimates of \code{phi} also increased by \emph{f} (this effect cannot be generally achieved if a non-trivial \code{resid.formula} with log link). Note that this is what \code{glm} does, but that some other widely used packages can behave differently (and obscurely). 
 }
  \item{processed}{
    A list of preprocessed arguments, for programming purposes only (as in \code{corrHLfit} code).
 }
}
\details{

 \bold{Fitting methods:}
 Many approximations for likelihood have been defined to fit mixed models (e.g. Noh and Lee (2007) for some overview), 
 and this function only considers a subset of them, but it adds a new complication in terms of REML methods. For example, 
 PQL as originally defined by Breslow and Clayton uses REML to estimate dispersion parameters, but this function allows one to use ML instead.
 Moreover, it allows some non-standard specification of the model formula that determines the conditional distribution used in REML.

 In the more general syntax for \code{HLmethod}, used as e.g. \code{HLmethod="RE(1,1)"}  
 the first '1' means that a first order Laplace approximation to the likelihood is used to estimate fixed effects 
 (a '0' would instead mean that the h likelihood is used as the objective function).
  The second  '1' means that a first order Laplace approximation to the likelihood or restricted likelihood 
  is used to estimate dispersion parameters, including the dv/d tau term specifically discussed by Lee & Nelder 2001, p. 997
  (a '0' would instead mean that these terms are ignored).
  
  It is possible to enforce the EQL approximation for estimation of dispersion parameter by adding a third index with value 0. 
  \code{"HL(0,1,0)"} is Lee & Nelder's (2001) method, i.e. \code{"EQL+"}.
  
  For a Gamma GLM with log link, ML and EQL results will differ in their phi estimates, and the EQL estimate will match that from the \code{\link{glm}} function. 
  

  \bold{Random effects} are constructed in several steps. first, a vector \bold{u} of independent and identically distributed (iid) random effects is drawn from some distribution;
 second, a transformation v=f(u) is applied to each element (this defines \bold{v} which elements are still iid); third, correlated random effects are obtained as \bold{Lv} 
  where \bold{L} is the \dQuote{square root} of a correlation matrix (this may be meaningful only for Gaussian random effects). 
  Finally, a matrix \bold{Z} (or sometimes \bold{ZA}, see \code{\link{Predictor}}) allows to specify how the correlated random effects
  affect the response values. In particular, \bold{Z} is the identity matrix if there is a single observation (response) for each location, but otherwise
  its elements \eqn{z_{ji}} are 1 for the \eqn{j}th observation in the \eqn{i}th location. 
  The design matrix for \bold{v} is then of the form \bold{ZL}. 

 The specification of the random effects u and v handles the following cases: 
 \describe{
 \item{Gaussian}{with zero mean, unit variance, and identity link;} 
 \item{Beta}{Beta-distributed random effects, where (\eqn{u ~ B(1/(2\lambda),1/(2\lambda))} 
    (mean=1/2, var\eqn{=\lambda/[4(1+\lambda)]}), with logit link \code{v=logit(u)})}
 \item{Gamma}{Gamma-distributed random effects, where \eqn{u ~ \Gamma(1/\lambda,\lambda)}, so u
   has mean 1 and variance \eqn{\lambda}. Both the log (\eqn{v=log(u)}) and identity (\eqn{v=u}) links are possible, though in the latter case the variance of \eqn{u} is constrained below 1. Gamma-distributed random effects with higher variance can be fitted using \code{inverse.Gamma(-1/mu)}. 
  } 
 \item{Inverse-Gamma}{inverse-Gamma distributed random effects, where \eqn{u ~ }inverse-Gamma\eqn{(1+1/\lambda,1/\lambda)}
    (mean=1, var=\eqn{\lambda/(1-\lambda))}. The allowed links are \code{v=log(u)} and \code{v=-1/u}. With the latter one, \code{v} 
    is Gamma-distributed with mean \eqn{\lambda} and variance \eqn{\lambda(\lambda+1)}.}
 }  

  \bold{The standard errors} reported may sometimes be misleading. For each set of parameters among \eqn{\beta}, \eqn{\lambda}, and \eqn{\phi} parameters these are computed assuming that the other parameters are known without error. This is why they are labelled \code{Cond. SE} (conditional standard error). This is most uninformative in the unusual case where \eqn{\lambda} and \eqn{\phi} are not separately estimable parameters. Further, the SEs for \eqn{\lambda} and \eqn{\phi} are rough approximations as discussed in particular by Smyth et al. (2001; \eqn{V_1} method).    

The extractor function \code{\link{get_any_IC}} can compute information criteria (\dQuote{AIC}) and effective degrees of freedom from the HLfit results. See the \code{get_any_IC} documentation for more details. 
}
\value{
An object of class \code{HLfit}, actually a list with many elements, several of which represent input arguments. 
Some elements may be undocumented. 

A few extractor functions are available (see \code{\link{extractors}}), 
and should be used as far as possible as they should be backward-compatible from version 1.4 onwards, while the structure of the return object may still evolve. The following information will be useful for extracting further elements of the object.

Elements describing the fit include \bold{descriptors of the fit}:

\item{eta}{Fitted values on the linear scale  (including the predicted random effects);}
\item{fv}{Fitted values (\eqn{\mu=}<inverse-link>(\eqn{\eta})) of the response variable (returned by the \code{fitted} function);}
\item{fixef}{The fixed effects coefficients, \eqn{\beta} (returned by the \code{fixef} function);}
\item{ranef}{The random effects \eqn{u} (returned by the \code{ranef} function);}
\item{v_h}{The random effects on the linear scale, \eqn{v};}
\item{phi}{The residual variance \eqn{\phi};}
\item{phi.object}{A possibly more complex object describing \eqn{\phi};}
\item{lambda}{The random effects (\eqn{u}) variance \eqn{\lambda};}
\item{lambda.object}{A possibly more complex object describing \eqn{\lambda};}
\item{corrPars}{Agglomerates information on correlation parameters, either fixed, or estimated by \code{HLfit} or \code{corrHLfit};}
\item{APHLs}{A list which elements are various likelihood components and information criteria. ikelihood components include conditional likelihood, the h-likelihood, and the two adjusted profile h-likelihoods:
   the (approximate) marginal likelihood \code{p_v} and the (approximate) restricted likelihood \code{p_bv} (the latter two available through the \code{logLik} function). Information criteria are described in Details and can be displayed by \code{\link{get_any_IC}};
}
\item{beta_cov}{Covariance matrix of \eqn{\beta} estimates.}

\bold{Information about the input} is contained in output elements named as \code{HLfit} or \code{corrHLfit} arguments (\code{data,family,resid.family,ranFix,prior.weights}), with the following notable exceptions or modifications:

\item{predictor}{The linear predictor, including the \code{formula} (possibly reformatted) and several attributes such as \code{ZALMatrix} which is the design matrix for random effects;}
\item{resid.predictor}{Analogous to \code{predictor}, for the residual variance;}
\item{rand.families}{corresponding to the \code{rand.family} input;}

\bold{Further miscellaneous diagnostics and descriptors of model structure:}

\item{leve_phi,lev_lambda}{Leverages;}
\item{ZAlist}{Components of the design matrix for random effects;}
\item{X.pv}{The design matrix for fixed effects;}
\item{fixef_terms,fixef_levels}{Futher information about fixed effect model;}
\item{weights}{(binomial data only) the binomial denominators;}
\item{y}{the response vector; for binomial data, the frequency response.}
\item{models}{Additional information on model structure for \eqn{\eta}, \eqn{\lambda} and \eqn{\phi};}
\item{HL}{A set of indices that characterize the approximations used for likelihood;}
\item{warnings}{A list of warnings for events that may have occurred during the fit.}

Finally, the object includes programming tools: \code{call, spaMM.version, get_w_h_coeffs, get_beta_w_cov, get_invColdoldList, get_logdispObject}.

}
\references{
Breslow, NE, Clayton, DG. (1993). Approximate Inference in Generalized Linear Mixed Models.
Journal of the American Statistical Association 88, 9-25.

Cox, D. R. and Donnelly C. A. (2011) Principles of Applied Statistics. Cambridge Univ. Press.

Ha, I. D., Lee, Y. and MacKenzie, G. (2007) Model selection for multi-component frailty models. Statistics in Medicine 26: 4790-4807.

Lee, Y., Nelder, J. A. (2001)  Hierarchical generalised linear models: A
synthesis of generalised linear models, random-effect models and structured
dispersions. Biometrika 88, 987-1006.

Lee, Y., Nelder, J. A. and Pawitan, Y. (2006). Generalised linear models with random effects: unified analysis via
h-likelihood. Chapman & Hall: London.

Noh, M., and Lee, Y. (2007). REML estimation for binary data in GLMMs, J.
Multivariate Anal. 98, 896-915.

Smyth GK, Huele AF, Verbyla AP (2001). Exact and approximate REML for heteroscedastic regression. Statistical Modelling 1, 161-175. 

Vaida, F., and Blanchard, S. 2005. Conditional Akaike information for mixed-effects models. Biometrika 92, 351-370.
}
\seealso{
\code{\link{HLCor}} for estimation with given spatial correlation parameters;
\code{\link{corrHLfit}} for joint estimation with spatial correlation parameters.
}

\examples{
data(wafers)
## Gamma GLMM with log link
HLfit(y ~X1+X2+X1*X3+X2*X3+I(X2^2)+(1|batch),family=Gamma(log),
          resid.model = ~ X3+I(X3^2) ,data=wafers)
%- : tested in update.Rd
## Gamma - inverseGamma HGLM with log link
HLfit(y ~X1+X2+X1*X3+X2*X3+I(X2^2)+(1|batch),family=Gamma(log),
          HLmethod="HL(1,1)",rand.family=inverse.Gamma(log),
          resid.model = ~ X3+I(X3^2) ,data=wafers)
}          
\keyword{ model }
