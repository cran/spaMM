\name{beta_resp}
\alias{beta_resp}
\title{
Beta-distribution family object
}
\description{
Returns a \code{family} object for beta-response models.
The precision parameter \code{prec} of this family is a positive value such that the variance of the response given its mean \eqn{\mu} is \eqn{\mu(1-\mu)/(1+}\code{prec}). \code{prec} is thus the precision parameter \eqn{\phi} of Ferrari & Cribari-Neto (2004) and of the \pkg{betareg} package (Cribari-Neto & Zeileis 2010).  

The family described by this object is not a GLM family, so it should not be used as a \code{\link[stats]{glm}} argument unless one additionally uses \code{method="llm.fit"} (see Examples for \code{\link{negbin1}} which is equivalent in this respect).
But it is still characterized by the same concepts as GLMs: a linear predictor, a link function, and the given distribution of residual variation. 

A \bold{fixed-effect} residual-dispersion model can be fitted, using the  \code{\link{resid.model}} argument, which is used to specify the form of the logarithm of the precision parameter (see Examples). Thus the variance of the response become \eqn{\mu(1-\mu)/(1+}\code{exp(<specified linear expression>)}).
}
\details{
Prior weights are meaningful for this family and handled as a factor of the precision parameter (as for GLM families) hence here not as a divisor of the variance (in contrast to GLM families): the variance of the response become \eqn{\mu(1-\mu)/(1+}\code{prec*<prior weights>}). However, this feature is experimental and may be removed in the future. The fitting function's \code{resid.model} argument may be preferred to obtain the same effect, by specifying an \code{offset(log(<prior weights>))} in its formula (given the log link used in that model). As usual in \pkg{spaMM}, the offset(.) argument should be a vector and any variable necessary for evaluating it should be in the \code{data}.
}
\usage{
beta_resp(prec, link = "logit")
}
\arguments{
  \item{prec}{
  Scalar (or left unspecified): precision parameter of the beta distribution. The \dQuote{outer-optimizing} functions, such as \code{fitme} (see \code{\link{inner-vs-outer}}) will estimate this parameter if called with \code{family=beta_prec()} (no given \code{prec} value). In other usages of this family, \code{prec} must be specified. 

}
  \item{link}{
  logit, probit, cloglog or cauchit link, specified by any of the available ways for GLM links (name, character string, one-element character vector, or object of class \code{link-glm} as returned by \code{\link[stats]{make.link}}). 
}
}
\value{
A list, formally of class \code{c("LLF", "family")}.  See \code{\link{LL-family}} for details about the structure and usage of such objects.
}
\references{
Cribari-Neto, F., & Zeileis, A. (2010). Beta Regression in R. Journal of Statistical Software, 34(2), 1-24. \doi{10.18637/jss.v034.i02}

Ferrari SLP, Cribari-Neto F (2004). “Beta Regression for Modelling Rates and Proportions.”
Journal of Applied Statistics, 31(7), 799-815.
}
\seealso{Further examples in \code{\link{LL-family}}.}
\examples{
  set.seed(123)
  beta_dat <- data.frame(y=runif(100),grp=sample(2,100,replace = TRUE), x_het=runif(100))
  
  fitme(y ~1+(1|grp), family=beta_resp(), data= beta_dat)
  ## same logL, halved 'prec' when prior weights=2 are used: 
  # fitme(y ~1+(1|grp), family=beta_resp(), data= beta_dat, prior.weights=rep(2,100))
  
  ## With model for residual dispersion:
  # fitme(y ~1+(1|grp), family=beta_resp(), data= beta_dat, resid.model= ~ x_het)
}
\keyword{models}
\keyword{regression}