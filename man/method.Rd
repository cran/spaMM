\name{method}
\alias{method}
\title{
Fitting methods (objective functions maximized)
}
\description{
 Many approximations for likelihood have been defined to fit mixed models (e.g. Noh and Lee (2007) for some overview), 
 and fitting functions in spaMM implement several of them, and some additional ones. In particular, 
 PQL as originally defined by Breslow and Clayton (1993) uses REML to estimate dispersion parameters, but spaMM allows one to use an ML variant of PQL. Moreover, it allows some non-standard specification of the model formula that determines the conditional distribution used in REML.
 
  EQL stands for the EQL method of Lee and Nelder (2001). The '+' version includes the d v/ d tau correction 
  described p. 997 of that paper, and the '-' version ignores it.
  PQL can be seen as the version of EQL- for GLMMs. It estimates fixed effects by maximizing h-likelihood and dispersion parameters by an approximation of REML, i.e. by maximization of an approximation of restricted likelihood.
  PQL/L is PQL without the leverage corrections that define REML estimation of random-effect parameters. Thus, it estimates dispersion parameters by an approximation of marginal likelihood.

 The \code{method} (or \code{HLmethod}) argument of fitting functions also accepts value of the form \code{"HL(<...>)"}, \code{"ML(<...>)"} and \code{"RE(<...>)"}, e.g. \code{method="RE(1,1)"}, which allow a more direct specification of the approximations used. 
 HL and RE are equivalent (both imply an REML correction).   
 The first '1' means that a first order Laplace approximation to the likelihood is used to estimate fixed effects 
 (a '0' would instead mean that the h likelihood is used as the objective function).
  The second  '1' means that a first order Laplace approximation to the likelihood or restricted likelihood 
  is used to estimate dispersion parameters, this approximation including the dv/d tau term specifically discussed by Lee & Nelder 2001, p. 997 (a '0' would instead mean that these terms are ignored).
  
  It is possible to enforce the EQL approximation for estimation of dispersion parameter (i.e., Lee and Nelder's (2001) method) by adding a third index with value 0. \code{"EQL+"} is thus \code{"HL(0,1,0)"}, while \code{"EQL-"} is \code{"HL(0,0,0)"}. \code{"PQL"} is EQL- for GLMMs. \code{"REML"} is \code{"HL(1,1)"}. \code{"ML"} is \code{"ML(1,1)"}. 
  
  Some of these distinctions make sense for \bold{GLMs}, and \code{\link[stats]{glm}} methods use approximations, which make a difference for Gamma GLMs. This means in particular that, (as stated in \code{stats::\link[stats]{logLik}}) the logLik of a Gamma GLM fit by \code{glm} differs from the exact likelihood. Further, the dispersion estimate returned by \code{summary.glm} differs from the one implied by \code{logLik}, because \code{summary.glm} uses Pearson residuals instead of deviance residuals. This may be confusing, and no \code{method} in spaMM tries to reproduce simultaneously these distinct features (however, \code{\link{spaMM_glm}} may do so). An \code{"ML(0,0,0)"} approximation of true ML provides the same log likelihood as \code{stats::logLik}, and the dispersion estimate returned by an \code{"HL(.,.,0)"} fit matches what can be computed from residual deviance and residual degrees of freedom of a glm fit, but this is not the estimate displayed by \code{summary.glm}. With a log link, the fixed effect estimates are unaffected by these distinctions.     
}
\references{
Breslow, NE, Clayton, DG. (1993). Approximate Inference in Generalized Linear Mixed Models.
Journal of the American Statistical Association 88, 9-25.

Lee, Y., Nelder, J. A. (2001)  Hierarchical generalised linear models: A
synthesis of generalised linear models, random-effect models and structured
dispersions. Biometrika 88, 987-1006.


Noh, M., and Lee, Y. (2007). REML estimation for binary data in GLMMs, J.
Multivariate Anal. 98, 896-915.
}
