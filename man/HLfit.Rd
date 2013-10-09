\name{HLfit}
\alias{HLfit}
\alias{inverse.Gamma}
\alias{AIC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fit mixed models with given correlation matrix}
\description{
  This fonction fits a class of GLMMs as well as some hierarchical generalized linear models (HGLM; Lee and Nelder 2001).
  HLfit fits both fixed effects parameters, and dispersion parameters i.e. the variance of the random effects 
  and the variance of the residual error. The linear predictor is of the 
  standard form \code{offset+ X beta + Z v}, where 
  X is the design matrix of fixed effects and Z is a design matrix of random effects. 
  The function also handles a linear predictor (with only fixed effects) for the residual variance. 
}
\usage{
HLfit(predictor, data, family = gaussian(), rand.family = gaussian(), 
      resid.predictor = ~1, REMLformula = NULL, verbose = c(F, F), 
      HLmethod = "HL(1,1)", control.HLfit = list(), init.HLfit = list(), 
      ranFix = list(), etaFix = list(), processed = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{predictor}{
  A \code{\link{Predictor}} including a \code{\link{formula}}, 
  or simply a \code{formula} if neither an offset nor a 
  spatial correlation matrix are required.
}
  \item{data}{
     A data frame containing the variables named in the model formula.  
}
  \item{family}{
   A \code{family} object describing the distribution of the response variable. 
   Possible values include the gaussian, poisson, binomial and Gamma families. Possible links are identity, log,
   inverse, logit, probit, and cloglog. Possible combinations of family and link are those allowed by the family objects
   (see \code{\link{family}}).
}
  \item{rand.family}{
  A \code{family} object describing the distribution of the random effect. Possible options are
  \code{gaussian()}, \code{Gamma(log)}, \code{Beta(logit)}, \code{inverse.Gamma(-1/mu)}, and \code{inverse.Gamma(log)}.
  For discussion of these alternatives see Lee and Nelder 2001 or Lee et al. 2006, p. 178-.
  Here the family gives the distribution of a random effect \eqn{u} 
  and the link gives \code{v} as function of \eqn{u} (see Details).      
}
  \item{resid.predictor}{
  A formula (without left-hand side) for the variance \code{phi} of the residual error. Currently can only contain fixed effects. 
}
  \item{REMLformula}{
  A model \code{formula} that allows the estimation of dispersion parameters, and 
  computation of restricted likelihood (\code{p_bv}) under a model different from the predictor \code{formula}.

  For example, if only random effects are included in \code{REMLformula}, an ML fit is performed and \code{p_bv} equals
  the marginal likelihood (or its approximation), \code{p_v}. This ML fit can be performed more simply by setting 
  \code{HLmethod="ML"} and leaving \code{REMLformula} at its default NULL value.
}
  \item{verbose}{
    A pair of booleans. The first controls messages about the results of the fit, the second controls messages about 
    the iterations.
  }
  \item{HLmethod}{
  Allowed values are \code{"REML"}, \code{"ML"}, \code{"REPQL"}, \code{"PQL/L"}, and expressions of the 
  form \code{"HL(<...>)"}, \code{"ML(<...>)"} and \code{"RE(<...>)"}. 
  The default behaviour is HL(1,1), which chooses between ML ('ML(1,1)') 
  and (possibly non-standard forms of) REML ('RE(1,1)'), according to \code{REMLformula}; 
  \bold{the joint default for these two arguments is REML}.

  The first '1' means that a first order Laplace approximation to the likelihood is used to estimate fixed effects.
  The second  '1' means that a first order Laplace approximation to the likelihood or restricted likelihood 
  is used to estimate dispersion parameters.
  'RE(0,1)' would mean that the h likelihood is used to estimate fixed effect parameters, as in PQL methods. 
  'REPQL' is Breslow and Clayton's PQL and could also be coded as 'RE(0,0)'. This is equivalent to 'RE(0,1)' except 
   that the latter involves an additional term in the approximation of the gradient of restricted likelihood. 
  'PQL/L' is likewise 'ML(0,0)' and equivalent to 'ML(0,1)' except that the latter involves the same additional term 
  but now for the approximation of the gradient of likelihood used for ML estimation of dispersion parameters. 
}
  \item{control.HLfit}{
  A list of parameters controlling (1) AIC computation; and (2) the fitting algorithms.

  \code{AIC=TRUE} provides the corrected AIC of Ha et al. 2007, 
  but the default is \code{FALSE} (think twice before applying model selection automatically). Note that the definition of AIC statistics 
  for mixed models is a current field of research.  

  Controls for the fitting algorithms should be ignored in routine use. They are
%
  \code{conv.threshold}, 
   a convergence threshold for the iterative algorithm, which controls whether 
   linear predictor terms (fixed effects and inferred random effects), and dispersion parameter estimates have converged. Defaults to 1e-05;
%
  \code{iter.mean.dispFix}, the number of iterations of the iterative algorithm for coefficients of the linear predictor,
       if no dispersion parameters are estimated by the iterative algorithm. Defaults to 40; 
%
  \code{iter.mean.dispVar}, the number of iterations of the iterative algorithm for coefficients of the linear predictor,
       if some dispersion parameter(s) is estimated by the iterative algorithm. Defaults to 20;  
%
  \code{max.iter}, the number of iterations of the iterative algorithm for joint estimation of dispersion parameters and
        of coefficients of the linear predictor. Defaults to 200.

}
  \item{init.HLfit}{
  A list of initial values for the iterative algorithm, with possible elements of the list are 
  \code{fixef} for fixed effect estimates (beta),  
  \code{v_h} for random effects (v, on linear scale),
  \code{lambda} for the variance of random effects (u), and 
  \code{phi} for the residual variance. 
  However, this argument can be ignored in routine use. 
}
  \item{ranFix}{
  A list of fixed values of random effect parameters, with possible elements \code{lambda}, and also \code{phi} for gaussian and Gamma HGLMs. 
  Inhibits the estimation of these parameters. 
 }
  \item{etaFix}{
   A list of fixed values of the coefficients of the linear predictor, with possible elements \code{beta} and \code{v_h} or \code{u_h} 
  (the latter two differ from each other for non-Gaussian random effects).
 }
  \item{processed}{
    A list of preprocessed arguments, for programming purposes only (as in \code{corrHLfit} code).
 }
}
\details{
 Many approximations for likelihood have been defined to fit mixed models (e.g. Noh and Lee (2007) for some overview), 
 and this function only considers a subset of them, but it adds a new complication in terms of REML methods. For example, 
 PQL as originally defined by Breslow and Clayton uses REML to estimate dispersion parameters, but this function allows one to use ML instead.
 Moreover, it allows some non-standard specification of the model formula that determines the conditional distribution used in REML.

  Random effects are constructed in several steps. first, independent and identically distributed (iid) random effects u are drawn from some distribution;
 second, a transformation v=f(u) is applied (v values are still iid); third, correlated random effects are obtained as \code{L \%*\% v} 
  where L is the \dQuote{square root} of a correlation matrix (this may be meaningful only for Gaussian random effects). 
  Finally, a matrix Z allows to specify how the correlated random effects
  affect the response values. In particular, Z is the identity matrix if there is a single observation (response) for each location, but otherwise
  its elements \eqn{z_{ji}} are 1 for each observation \eqn{j} in a given location \eqn{i}. 
  The design matrix for the \code{v} is then of the form \code{Z \%*\% L}. 

 The specification of the random effects u and v handles the following cases: 
 \describe{
 \item{Gaussian}{with zero mean, unit variance, and identity link;} 
 \item{Beta}{Beta-distributed random effects, where (\eqn{u ~ B(1/(2\lambda),1/(2\lambda))} 
    (mean=1/2, var\eqn{=\lambda/[4(1+\lambda)]}), with logit link \code{v=logit(u)})}
 \item{Gamma}{Gamma-distributed random effects, where \eqn{u ~ \Gamma(\lambda,1/\lambda)}, so u
   has mean 1 and variance \eqn{\lambda}. With a log link (\eqn{v=log(u ~ \Gamma)}), v has mean zero.} 
 \item{Inverse-Gamma}{inverse-Gamma distributed random effects, where \eqn{u ~ }inverse-Gamma\eqn{(1+1/\lambda,1/\lambda)}
    (mean=1, var=\eqn{\lambda/(1-\lambda))}, with \code{v=log(u)} or \code{v=-1/u}.}
 }

}
\value{
An object of class \code{HLfit}, actually a list with many elements, many of which represent input arguments. The \code{summary.HLfit} function 
normally formats this nicely, so that the following information is useful mainly if one wants to extract any of the elements of the object. 
These elements include, first, variables describing the fit:  

\item{APHLs}{A list with usually four elements, the conditional likelihood, 
   the h-likelihood, and the two adjusted profile h-likelihoods:
   the (approximate) marginal likelihood \code{p_v} and the (approximate) restricted likelihood \code{p_bv}.
}
\item{fixef}{The fixed effects coefficients, \eqn{\beta}}
\item{fixef_se}{Estimates of standard errors for \eqn{\beta} estimates}
\item{ranef}{The random effects \eqn{u}}
\item{v_h}{The random effects on the linear scale, \eqn{v}}
\item{phi}{The residual variance \eqn{\phi}}
\item{phi.object}{A possibly more complex object describing \eqn{\phi}}
\item{lambda}{The random effects (\eqn{u}) variance \eqn{\lambda}}
\item{lambda.object}{A possibly more complex object describing \eqn{\lambda}}
\item{eta}{Fitted values on the linear scale  (including the predicted random effects)}
\item{fv}{Fitted values (\eqn{\mu=}<inverse-link>(\eqn{\eta})) of the response variable }

It may also be worth checking the following elements 

\item{HL}{A set of indices that characterize the approximations used for likelihood.}
\item{warnings}{A list of warnings of events that may have occurred during the fit.}

Additional information about the fit is contained in

\item{ZLmatrix}{The design matrix for random effects (see Details).}

and in (possibly reformatted) input arguments

\item{predictor}{The input predictor.}
\item{data}{The input data.}
\item{family}{a \code{\link{family}} object corresponding to the \code{family} input}
\item{rand.family }{a \code{\link{family}} object corresponding to the \code{rand.family} input}
\item{y}{the response vector; for binomial data, the frequency response}
\item{formulaDisp}{Model formula for dispersion response}
\item{RespLink_disp}{Link for dispersion response}
\item{X}{The design matrix for fixed effects}
\item{ranFix}{A \code{ranFix} input}
\item{corrPars}{Additional information on correlation parameters, not necessarily used by HLfit itself but in upper calling functions such as HLCor or corrHLfit}
\item{models}{Additional information on model structure for \eqn{\eta}, \eqn{\lambda} and \eqn{\phi}}
\item{weights}{(binomial data only) the binomial denominators}
}
\references{
Ha, I. D., Lee, Y. and MacKenzie, G. (2007) Model selection for multi-component frailty models. Statistics in Medicine 26: 4790-4807.

Lee, Y., Nelder, J. A. (2001)  Hierarchical generalised linear models: A
synthesis of generalised linear models, random-effect models and structured
dispersions. Biometrika 88, 987-1006.

Lee, Y., Nelder, J. A. and Pawitan, Y. (2006). Generalised linear models with random effects: unified analysis via
h-likelihood. Chapman & Hall: London.

Noh, M., and Lee, Y. (2007). REML estimation for binary data in GLMMs, J.
Multivariate Anal. 98, 896-915.
}
\seealso{
\code{\link{corrHLfit}} for joint estimation with correlation parameters.
}

\examples{
data(wafers)
## Gamma GLMM with log link (also called Gamma-Lognormal model)
HLfit(predictor= y ~X1+X2+X1*X3+X2*X3+I(X2^2)+(1|batch),family=Gamma(log),
          resid.predictor = ~ X3+I(X3^2) ,data=wafers)
## Gamma - inverseGamma HGLM with log link
HLfit(predictor= y ~X1+X2+X1*X3+X2*X3+I(X2^2)+(1|batch),family=Gamma(log),
          HLmethod="HL(1,1)",rand.family=inverse.Gamma(log),
          resid.predictor = ~ X3+I(X3^2) ,data=wafers)
}
\keyword{ model }
