\name{COMPoisson}
\alias{COMPoisson}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Conway-Maxwell-Poisson (COM-Poisson) GLM family
}
\description{
The COM-Poisson family is a generalization of the Poisson family which can describe over-dispersed as well as under-dispersed count data. It is indexed by a parameter \eqn{\nu_{CMP}} that quantifies such dispersion. The COM-Poisson family is here implemented as a \code{\link[stats]{family}} object, so that it can be fitted by \code{\link[stats]{glm}}, and further used to model conditional responses in mixed models fitted by this package's functions (see Examples). \eqn{\nu_{CMP}} is distinct from the dispersion parameter \eqn{\nu=1/\phi} considered elsewhere in this package and in the GLM literature, as \eqn{\nu} affects in a more specific way the log-likelihood. The \dQuote{canonical link} \eqn{\theta(\mu)} between the canonical GLM parameter \eqn{\theta} and the expectation \eqn{\mu} of the response does not have a known expression in terms of elementary functions. The link inverse is \eqn{\mu(\theta)=\sum_(i=0)^\infty \lambda^i / (i!)^\nu} for \eqn{\lambda=e^\theta} (hence the link is here nicknamed \code{"loglambda"}). 
}
\usage{
COMPoisson(nu =  stop("COMPoisson's 'nu' must be specified"), 
           link = "loglambda")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link}{
GLM link function. Cannot be modified.
}
  \item{nu}{
  Under-dispersion parameter \eqn{\nu_{CMP}}. The \code{fitme} and \code{corrHLfit} functions called with \code{family=COMPoisson()} (no given \eqn{\nu_{CMP}} value) will estimate \eqn{\nu_{CMP}}. In other usage, \eqn{\nu_{CMP}} must be specified. \code{COMPoisson(nu=1)} is the Poisson family. For \eqn{\nu>1}, the distribution is under-dispersed. The limit as \eqn{\nu -> 0} is the Bernoulli distribution with expectation \eqn{\lambda/(1+\lambda)}. \code{COMPoisson(nu=0)} would describe a geometric distribution with parameter \eqn{\lambda}. However the code may fail to fit distributions with \eqn{\nu} approaching 0 because it evaluates various sums that have to be truncated and cannot be easily approximated. However it is possible to fit the case \eqn{\nu=0}, thanks to an ad-hoc algorithm for this case.     
}
}
\value{
A family object.
}
\details{
Fitting will be very slow when \eqn{\nu_{CMP}} approaches zero, i.e. for strong overdispersion of (conditional) response.
It may be faster if the value of the control parameter \code{COMP_maxn} is reduced by using \code{spaMM.options()}, but this will reduce the accuracy of the results for high overdispersion. This trade-off occurs because the link inverse function, as shown in Description, involves an infinite summation. This sum can be easily approximated for large \eqn{\nu_{CMP}} but not when \eqn{\nu_{CMP}} approaches zero. By default, \code{spaMM} truncates the sum at \code{spaMM.getOption(COMP_maxn)=} 10000 terms, which should be more than enough for underdispersed response. }
\references{
  G. Shmueli, T. P. Minka, J. B. Kadane, S. Borle and P. Boatwright (2005) A useful distribution for fitting discrete data: revival of the Conway-Maxwell-Poisson distribution. Appl. Statist. 54: 127-142.

	Sellers KF, Shmueli G (2010) A Flexible Regression Model for Count Data. Ann. Appl. Stat. 4: 943–961
}
\examples{
if (spaMM.getOption("example_maxtime")>7) {          
  # Fitting COMPoisson model with estimated nu parameter:
  data(freight) ## example from Sellers & Shmueli, Ann. Appl. Stat. 4: 943–961 (2010)
  fitme(broken ~ transfers, data=freight, family = COMPoisson())
  # GLMM with under-dispersed conditional response
  HLfit(broken ~ transfers+(1|id), data=freight, family = COMPoisson(nu=10),HLmethod="ML")
}
\dontrun{
  data(freight)
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