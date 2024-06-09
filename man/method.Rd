\name{method}
\alias{method}
\alias{obsInfo}
\title{
Fitting methods (objective functions maximized)
}
\description{
The \code{method} argument of the fitting functions, with possible values \code{"ML"}, \code{"REML"},\code{"PQL"}, \code{"PQL/L"}, and so on, controls whether restricted likelihood techniques are used to estimate residual variance and random-effect parameters, and the way likelihood functions are approximated.

By default, Laplace approximations are used, as selected by \code{"ML"} and \code{"REML"} methods. The Laplace approximation to (log-)marginal likelihood can be expressed in terms of the joint log-likelihood of the data and the random effects (or the \emph{h}-likelihood in works of Lee and Nelder). The Laplace approximation is the joint likelihood minus half the log-determinant of the matrix of second derivatives (Hessian matrix) of the negative joint likelihood with respect to the random effects (observed information matrix). The Laplace approximation to restricted likelihood (for REML estimation) is similarly defined from the Hessian matrix with respect to random effects \bold{and} fixed effects (for the adventurous, \pkg{spaMM} allows some non-standard specification of the fixed effects included in the definition of he Hessian).
  
Various additional approximations have been considered. Penalized quasi-likelihood (PQL), as originally defined for GLMMs by Breslow and Clayton (1993), uses a Laplace approximation of restricted likelihood to estimate dispersion parameters, and estimates fixed effects by maximizing the joint likelihood (h-likelihood). Although PQL has been criticized as an approximation of likelihood (and actual implementations may diverge from the original version), it has some interesting inferential properties. \pkg{spaMM} allows one to use an ML variant of PQL, named PQL/L.  

Further approximations defined by Lee, Nelder and collaborators (e.g., Noh and Lee, 2007, for some overview) may mostly be seen as laying between PQL and the full Laplace method in terms of approximation of the likelihood, and as extending them to models with non-gaussian random effects (\dQuote{HGLMs}). 
In practice the ML, REML, PQL and PQL/L methods should cover most (all?) needs for GLMMs, and EQL extends the PQL concept to HGLMs. \code{method="EQL+"} stands for the EQL method of Lee and Nelder (2001). The '+' signals that it includes the d v/ d tau correction described p. 997 of that paper, while \code{method="EQL-"} ignores it. \code{"PQL"} is equivalent to \code{EQL-} for GLMMs. \code{"PQL/L"} is PQL without the leverage corrections that characterize REML estimation of random-effect parameters. 

 \bold{\pkg{spaMM} uses the observed information matrix by default since version 4.0.0}. By contrast, in Laplace approximations of likelihood described in the work of Lee & Nelder, i.e. for mixed-effect models with GLM response families, the information matrix is written in terms of the GLM weights (e.g., Lee & Nelder 2001, p.1004), and is thus effectively the expected information matrix, which differs from the observed information matrix in the case of GLM families with non-canonical link (McCullagh & Nelder 1989, p.42). Therefore, the likelihood approximation based on the expected information matrix differs from the one based on the observed information matrix in the same conditions. 

For non-GLM response families (currently, the \code{\link{negbin1}}, \code{\link{beta_resp}} and \code{\link{betabin}}), only observed information is available (expected information would at best be quite difficult to evaluate, with no benefits). For GLM response families, use of expected information matrix can be required at a global level by setting \code{spaMM.options(obsInfo=FALSE)} or in a specific fit by adding \code{"exp"} as a second specifier in the method (e.g., \code{method=c("ML","exp")}). This can be distinctly useful (in terms of speed) for fitting models with \code{Gamma(log)} response family. Conversely, the \code{"obs"} specifier will enforce use of observed information matrix when the alternative is set at a global level.
}

\details{
The \code{method} (or \code{HLmethod}) argument of fitting functions also accepts values of the form \code{"HL(<...>)"}, \code{"ML(<...>)"} and \code{"RE(<...>)"}, e.g. \code{method="RE(1,1)"}, which allow one to experiment with further combinations of approximations. HL and RE are equivalent (both imply an REML correction). The first '1' means that a Laplace approximation to the likelihood is used to estimate fixed effects 
 (a '0' would instead mean that the h likelihood is used as the objective function). The second  '1' means that a Laplace approximation to the likelihood or restricted likelihood is used to estimate dispersion parameters, this approximation including the dv/d tau term specifically discussed by Lee & Nelder 2001, p. 997 (a '0' would instead mean that these terms are ignored). It is possible to enforce the EQL approximation for estimation of dispersion parameter (i.e., Lee and Nelder's (2001) method) by adding a third index with value 0. \code{"EQL+"} is thus \code{"HL(0,1,0)"}, while \code{"EQL-"} is \code{"HL(0,0,0)"}. \code{"PQL"} is \code{EQL-} for GLMMs. \code{"REML"} is \code{"HL(1,1)"}. \code{"ML"} is \code{"ML(1,1)"}. 
  
Some of these distinctions make sense for \bold{GLMs}, and may help in understanding idiosyncrasies of \code{stats::\link[stats]{glm}} for Gamma GLMs. In particular (as stated in the \code{stats::\link[stats]{logLik}} documentation) the logLik of a Gamma GLM fit by \code{glm} differs from the exact likelihood. An \code{"ML(0,0,0)"} approximation of true ML provides the same log likelihood as \code{stats::logLik}. Further, the dispersion estimate returned by \code{summary.glm} differs from the one implied by \code{logLik}, because \code{summary.glm} uses Pearson residuals instead of deviance residuals. This may be confusing, and no \code{method} in \pkg{spaMM} tries to reproduce simultaneously these distinct features (however, \code{\link{spaMM_glm}} may do so). The dispersion estimate returned by an \code{"HL(.,.,0)"} fit matches what can be computed from residual deviance and residual degrees of freedom of a \code{glm} fit, but this is not the estimate displayed by \code{summary.glm}. The fixed effect estimates are not affected by these tinkerings.     
}

\references{
Breslow, NE, Clayton, DG. (1993). Approximate Inference in Generalized Linear Mixed Models.
Journal of the American Statistical Association 88, 9-25.

Lee, Y., Nelder, J. A. (2001)  Hierarchical generalised linear models: A
synthesis of generalised linear models, random-effect models and structured
dispersions. Biometrika 88, 987-1006.

McCullagh, P. and Nelder, J.A. (1989) Generalized Linear Models, 2nd edition. London: Chapman & Hall.

Noh, M., and Lee, Y. (2007). REML estimation for binary data in GLMMs, J.
Multivariate Anal. 98, 896-915.
}
