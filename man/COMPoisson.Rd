\name{COMPoisson}
\alias{COMPoisson}
\alias{geometric}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Conway-Maxwell-Poisson (COM-Poisson) GLM family
}
\description{
The COM-Poisson family is a generalization of the Poisson family which can describe over-dispersed as well as under-dispersed count data. It is indexed by a parameter \code{nu} that quantifies such dispersion. It includes the Poisson, geometric and Bernoulli as special (or limit) cases (see Details). The COM-Poisson family is here implemented as a \code{\link[stats]{family}} object, so that it can be fitted by \code{\link[stats]{glm}}, and further used to model conditional responses in mixed models fitted by this package's functions (see Examples). \code{nu} is distinct from the dispersion parameter \eqn{\nu=1/\phi} considered elsewhere in this package and in the GLM literature, as \eqn{\nu} affects in a more specific way the log-likelihood. The \dQuote{canonical link} \eqn{\theta(\mu)} between the canonical GLM parameter \eqn{\theta} and the expectation \eqn{\mu} of the response does not have a known expression in terms of elementary functions. The link inverse is \eqn{\mu(\theta)=\sum_(i=0)^\infty \lambda^i / (i!)^\nu} for \eqn{\lambda=e^\theta} (hence the link is here nicknamed \code{"loglambda"}). 
}
\usage{
COMPoisson(nu =  stop("COMPoisson's 'nu' must be specified"), 
           link = "loglambda")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link}{GLM link function. Cannot be modified.}
  \item{nu}{
  Under-dispersion parameter. The \code{fitme} and \code{corrHLfit} functions called with \code{family=COMPoisson()} (no given \code{nu} value) will estimate this parameter. In other usage of this family, \code{nu} must be specified. \code{COMPoisson(nu=1)} is the Poisson family. }
}
\value{
A family object.
}
\details{
For \code{nu}>1, the distribution is under-dispersed. The limit as \code{nu}\eqn{\to \infty} is the Bernoulli distribution with expectation \eqn{\lambda/(1+\lambda)}. 

The link inverse function, as shown in Description, involves an infinite summation. In this summation and related computations for the COMPoisson model, the sum can be easily approximated by a finite sum for large \code{nu} but not when \code{nu} approaches zero. For this reason, the code may fail to fit distributions with \code{nu} approaching 0 (strong residual over-dispersion). The case \code{nu=0} itself is the geometric distribution with parameter \eqn{\lambda} and is fitted by an ad hoc algorithm devoid of such problems. Otherwise, \code{spaMM} truncates the sum, and uses numerical integrals to approximate missing terms (which slows down the fitting operation). In addition, it applies an ad hoc continuity correction to ensure continuity of the result in \code{nu=1} (Poisson case). These corrections affect numerical results for the case of residual overdispersion but are negligible for the case of residual underdispersion. Alternatively, \code{spaMM} uses Gaunt et al.'s approximations when the condition defined in \code{spaMM.getOption("CMP_asympto_cond")} is satisfied. 

The name \code{COMP_nu} should be used to set values of \code{nu} in control arguments of the fitting functions (e.g., \code{fitme(.,init=list(COMP_nu=1))}). 
}
\references{
  Gaunt, Robert E. and Iyengar, Satish and Olde Daalhuis, Adri B. and Simsek, Burcin. An asymptotic expansion for the normalizing constant of the Conway--Maxwell--Poisson distribution. Ann Inst Stat Math (2017) \doi{10.1007/s10463-017-0629-6}.
  
  G. Shmueli, T. P. Minka, J. B. Kadane, S. Borle and P. Boatwright (2005) A useful distribution for fitting discrete data: revival of the Conway-Maxwell-Poisson distribution. Appl. Statist. 54: 127-142.

	Sellers KF, Shmueli G (2010) A Flexible Regression Model for Count Data. Ann. Appl. Stat. 4: 943–961
}
\examples{
# Fitting COMPoisson model with estimated nu parameter:
data("freight") ## example from Sellers & Shmueli, Ann. Appl. Stat. 4: 943–961 (2010)
fitme(broken ~ transfers, data=freight, family = COMPoisson())
# GLMM with under-dispersed conditional response
HLfit(broken ~ transfers+(1|id), data=freight, family = COMPoisson(nu=10),HLmethod="ML")

\dontrun{
  data("freight")
  # Equivalence of poisson() and COMPoisson(nu=1):
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