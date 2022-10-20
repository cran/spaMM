\name{COMPoisson}
\alias{COMPoisson}
\alias{geometric}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Conway-Maxwell-Poisson (COM-Poisson) GLM family
}
\description{
The COM-Poisson family is a generalization of the Poisson family which can describe over-dispersed as well as under-dispersed count data. It is indexed by a parameter \code{nu} that quantifies such dispersion. For \code{nu}>1, the distribution is under-dispersed relative to the Poisson distribution with same mean. It includes the Poisson, geometric and Bernoulli as special (or limit) cases (see Details). The COM-Poisson family is here implemented as a \code{\link[stats]{family}} object, so that it can be fitted by \code{\link[stats]{glm}}, and further used to model conditional responses in mixed models fitted by this package's functions (see Examples). \code{nu} is distinct from the dispersion parameter \eqn{\nu=1/\phi} considered elsewhere in this package and in the GLM literature, as \eqn{\nu} affects in a more specific way the log-likelihood. 

Several links are now allowed for this family, corresponding to different versions of the COMPoisson described in the literature (e.g., Sellers & Shmueli 2010; Huang 2017).
}
\usage{
COMPoisson(nu =  stop("COMPoisson's 'nu' must be specified"), 
           link = "loglambda")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link}{GLM link function. The default is the canonical link \code{"loglambda"} (see Details), but other links are allowed (currently log, sqrt or identity links as commonly handled for the Poisson family).}
  \item{nu}{
  Under-dispersion parameter. The \code{fitme} and \code{corrHLfit} functions called with \code{family=COMPoisson()} (no given \code{nu} value) will estimate this parameter. In other usage of this family, \code{nu} must be specified. \code{COMPoisson(nu=1)} is the Poisson family. }
}
\value{
A family object.
}
\details{
The \eqn{i}th term of the distribution can be written \eqn{q_i/Z} where \eqn{q_i=\lambda^i / (i!)^\nu} and \eqn{Z=\sum_{(i=0)}^\infty q_i}, for \eqn{\lambda=\lambda(\mu)} implied by its inverse relationship, the expectation formula \eqn{\mu=\mu(\lambda)=\sum_{(i=0)}^\infty i q_i(\lambda)/Z}. The case \code{nu=0} is the geometric distribution with parameter \eqn{\lambda}; \code{nu=1} is the Poisson distribution with mean \eqn{\lambda}; and the limit as \code{nu} -> \eqn{\infty} is the Bernoulli distribution with expectation \eqn{\lambda/(1+\lambda)}. 

From this definition, this is an exponential family model with canonical parameters \eqn{log(\lambda)} and \eqn{\nu}. When the linear predictor \eqn{\eta} specifies \eqn{log(\lambda(\mu))}, the canonical link is used (e.g., Sellers & Shmueli 2010). It is here nicknamed \code{"loglambda"} and does not have a known expression in terms of elementary functions. To obtain \eqn{\mu} as the link inverse of the linear predictor \eqn{\eta}, one then first computes \eqn{\lambda=e^\eta} and then \eqn{\mu(\lambda)} by the  expectation formula. For other links (Huang 2017), one directly computes \eqn{\mu} by the link inverse (e.g., \eqn{\mu=e^\eta} for link \code{"log"}), and then one may solve for \eqn{\lambda= \lambda(\mu)} to obtain other features of the distribution.    

The relationships between \eqn{\lambda} and \eqn{\mu} or other moments of the distribution involve infinite summations. These sums can be easily approximated by a finite number of terms for large \code{nu} but not when \code{nu} approaches zero. For this reason, the code may fail to fit distributions with \code{nu} approaching 0 (strong residual over-dispersion). The case \code{nu=0} (the geometric distribution) is fitted by an ad hoc algorithm devoid of such problems. Otherwise, \code{spaMM} truncates the sum, and uses numerical integrals to approximate missing terms (which slows down the fitting operation). In addition, it applies an ad hoc continuity correction to ensure continuity of the result in \code{nu=1} (Poisson case). These corrections affect numerical results for the case of residual overdispersion but are negligible for the case of residual underdispersion. Alternatively, \code{spaMM} uses Gaunt et al.'s (2017) approximations when the condition defined in \code{spaMM.getOption("CMP_asympto_cond")} is satisfied. All approximations reduces the accuracy of computations, in a way that can impede the extended Levenberg-Marquardt algorithm sometimes needed by spaMM.

The name \code{COMP_nu} should be used to set initial values or bounds on \code{nu} in control arguments of the fitting functions (e.g., \code{fitme(.,init=list(COMP_nu=1))}). Fixed values should be set by the family argument (\code{COMPoisson(nu=.)}). 
}
\references{
  Gaunt, Robert E. and Iyengar, Satish and Olde Daalhuis, Adri B. and Simsek, Burcin. (2017) An asymptotic expansion for the normalizing constant of the Conway--Maxwell--Poisson distribution. Ann Inst Stat Math \doi{10.1007/s10463-017-0629-6}.
  
  Huang, Alan (2017) Mean-parametrized Conway-Maxwell-Poisson regression models for dispersed counts. Stat. Modelling \doi{10.1177/1471082X17697749} 
  
  G. Shmueli, T. P. Minka, J. B. Kadane, S. Borle and P. Boatwright (2005) A useful distribution for fitting discrete data: revival of the Conway-Maxwell-Poisson distribution. Appl. Statist. 54: 127-142.

	Sellers KF, Shmueli G (2010) A Flexible Regression Model for Count Data. Ann. Appl. Stat. 4: 943–961
}
\examples{
if (spaMM.getOption("example_maxtime")>0.9) {
  # Fitting COMPoisson model with estimated nu parameter:
  #
  data("freight") ## example from Sellers & Shmueli, Ann. Appl. Stat. 4: 943–961 (2010)
  fitme(broken ~ transfers, data=freight, family = COMPoisson())
  fitme(broken ~ transfers, data=freight, family = COMPoisson(link="log"))

  # glm(), HLCor() and HLfit() handle spaMM::COMPoisson() with fixed overdispersion:
  #
  glm(broken ~ transfers, data=freight, family = COMPoisson(nu=10))
  HLfit(broken ~ transfers+(1|id), data=freight, family = COMPoisson(nu=10),method="ML")
  
  # Equivalence of poisson() and COMPoisson(nu=1):
  #
  COMPglm <- glm(broken ~ transfers, data=freight, family = poisson())
  coef(COMPglm)
  logLik(COMPglm)
  COMPglm <- glm(broken ~ transfers, data=freight, family = COMPoisson(nu=1))
  coef(COMPglm)
  logLik(COMPglm)
  HLfit(broken ~ transfers, data=freight, family = COMPoisson(nu=1))

}
}
\keyword{models}
\keyword{regression}