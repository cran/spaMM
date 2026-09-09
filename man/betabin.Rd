\name{betabin}
\alias{betabin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Beta-binomial family object
}
\description{
Returns a \code{family} object for beta-binomial models.
The precision parameter \code{prec} of this family is a positive value such that the variance of the beta-distributed latent variable given its mean \eqn{\mu} is \eqn{\mu(1-\mu)/(1+}\code{prec}). \code{prec} is thus the same precision parameter as for the beta family (see \code{\link{beta_resp}}. The variance of the beta-binomial sample of size \eqn{n} is response is \eqn{\mu(1-\mu)n(n+}\code{prec}\eqn{)/(1+}\code{prec}). 

The family described by this object is not a GLM family, so it should not be used as a \code{\link[stats]{glm}} argument unless one additionally uses \code{method="llm.fit"} (see Examples).
But it is still characterized by the same concepts as GLMs: a linear predictor, a link function, and the given distribution of residual variation. 

A \bold{fixed-effect} residual-dispersion model can be fitted, using the  \code{\link{resid.model}} argument, which is used to specify the form of the logarithm of the precision parameter (see Examples). Thus the variance of the latent beta-distributed variable becomes \eqn{\mu(1-\mu)/(1+}\code{exp(<specified linear expression>)}).
}
\usage{
betabin(prec, link = "logit")
}
\details{
Prior weights are meaningful for this family and handled as a factor of the precision parameter of the latent beta-distributed variable: the variance of the latent variable become \eqn{\mu(1-\mu)/(1+}\code{prec*<prior weights>}). However, this feature is experimental and may be removed in the future. The fitting function's \code{resid.model} argument may be preferred to obtain the same effect, by specifying an \code{offset(log(<prior weights>))} in its formula (given the log link used in that model). As usual in \pkg{spaMM}, the offset(.) argument should be a vector and any variable necessary for evaluating it should be in the \code{data}.
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prec}{
  Scalar (or left unspecified): precision parameter of the beta distribution. The \dQuote{outer-optimizing} functions, such as \code{fitme} (see \code{\link{inner-vs-outer}}) will estimate this parameter if called with \code{family=betabin()} (no given \code{prec} value). In other usages of this family, \code{prec} must be specified. 
}
  \item{link}{
  logit, probit, cloglog or cauchit link, specified by any of the available ways for GLM links (name, character string, one-element character vector, or object of class \code{link-glm} as returned by \code{\link[stats]{make.link}}). 
}
}
\value{
A list, formally of class \code{c("LLF", "family")}.  See \code{\link{LL-family}} for details about the structure and usage of such objects.
}
\examples{
  data(clinics)
  fitme(cbind(npos,nneg)~1+(1|clinic), family=betabin(), data=clinics)
  
  # Compatibility with glm(., method="llm.fit"):
  glm(cbind(npos,nneg)~1, family=betabin(2), data=clinics, 
      method="llm.fit")
  # (trying to use default method=stats::glm.fit in this case will fail).
}
\keyword{models}
\keyword{regression}