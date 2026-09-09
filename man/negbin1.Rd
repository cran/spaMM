\name{negbin1}
\alias{negbin1}
\title{
Alternative negative-binomial family  
}
\description{
Returns a \code{family} object suitable as a \code{fitme} argument for fitting negative-binomial models with variance linearly (affinely) related to the mean \eqn{\mu}: variance=\eqn{\mu+\mu}/shape. This distribution can be represented as a Poisson-Gamma mixture, where the conditional Poisson mean is \eqn{\mu} times a Gamma random variable with mean 1 and variance 1/(\code{shape*}\eqn{\mu}) as produced by \code{rgamma(., shape=sh,scale=1/sh)} where \code{sh=shape*}\eqn{\mu}, meaning that the family \code{shape} parameter controls, but differs from, the gamma \code{shape} parameter. 

The family described by this object is not a GLM family, so it should not be used as a \code{\link[stats]{glm}} argument unless one additionally uses \code{method="llm.fit"} (see Examples).
But it is still characterized by the same concepts as GLMs: a linear predictor, a link function, and the given distribution of residual variation. 

The zero-truncated variant of this family is also handled, and is never suitable as argument for \code{glm}.

A \bold{fixed-effect} residual-dispersion model can be fitted, using the \code{\link{resid.model}} argument, which is used to specify the form of the logarithm of the shape parameter. Thus the variance of the response become \eqn{\mu+\mu}/\code{exp(<specified linear expression>)}.
}
\usage{
negbin1(shape, link = "log", trunc = -1L)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{shape}{
  Shape parameter. The \dQuote{outer-optimizing} functions, such as \code{fitme} (see \code{\link{inner-vs-outer}}) will estimate this parameter if called with \code{family=negbin1()} (no given \code{shape} value). In other usages of this family, \code{shape} must be specified. 
}
  \item{link}{
  log, sqrt or identity link, specified by any of the available ways for GLM links (name, character string, one-element character vector, or object of class \code{link-glm} as returned by \code{\link[stats]{make.link}}). 
}
  \item{trunc}{
  Either \code{0L} for zero-truncated distribution, or \code{-1L} for default untruncated distribution.
  }
}
\details{
The name \code{NB_shape} should be used to set values of shape in optimization control arguments of the fitting functions (e.g., \code{fitme(.,init=list(NB_shape=1))}); but fixed values are set by the \code{shape} argument.

% For approx as negbin2 with an heteroscedastic Gamma random-effect, specified e.g. as \code{(weights-1|.)} where the weights need to be updated iteratively as function of predicted response, cf file \code{test-negbin1.R} in the \code{/test} directory provides one example. Other mean-variance relationship can be handled in the same way.
}
\value{
A list, formally of class \code{c("LLF", "family")}. See \code{\link{LL-family}} for details about the structure and usage of such objects.
}
\seealso{Examples in \code{\link{LL-family}}. \code{\link{resid.model}} for an example with a residual-dispersion model.}
\examples{
data(scotlip)
## Shape estimated:
fitme(I(1+cases)~I(prop.ag/10)+offset(log(expec)),family=negbin1(), data=scotlip)

## Shape fixed:
fitme(I(1+cases)~I(prop.ag/10)+offset(log(expec)),family=negbin1(1), data=scotlip)

# => consistent with:
glm(I(1+cases)~I(prop.ag/10)+offset(log(expec)),family=negbin1(1), 
    data=scotlip, method="llm.fit")

# !! but not with !!
# glm(I(1+cases)~I(prop.ag/10)+offset(log(expec)),family=negbin1(1), 
#     data=scotlip) # WRONG !!
}
\keyword{models}
\keyword{regression}