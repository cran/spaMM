\name{negbin1}
\alias{negbin1}
\title{
Alternative negative-binomial family  
}
\description{
Returns a \code{family} object suitable as a \code{fitme} argument for fitting negative-binomial models with variance linearly (affinely) related to the mean \eqn{\mu}: variance=\eqn{\mu+\mu}/shape, with known or unknown underlying Gamma shape parameter. 
The model described by such a family is characterized by a linear predictor, a link function, and such a negative-binomial model for the residual variation. The zero-truncated variant of this family is also handled.

A \bold{fixed-effect} residual-dispersion model can be fitted, using the \code{\link{resid.model}} argument, which is used to specify the form of the logarithm of the shape parameter. Thus the variance of the response become \eqn{\mu+\mu}/\code{exp(<specified linear expression>)}.
}
\usage{
negbin1(shape = stop("negbin1's 'shape' must be specified"), link = "log", trunc = -1L)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{shape}{
  Parameter controlling the mean-variance relationship of the negative binomial distribution. Given that the distribution can be represented as a Poisson-Gamma mixture, where the conditional Poisson mean is \eqn{\mu} times a Gamma random variable with mean 1 and (distinct) shape \code{sh}  (as produced by \code{rgamma(., shape=sh,scale=1/sh)}), the latter shape is then determined as \code{sh=shape*}\eqn{\mu}. 
}
  \item{link}{
  log, sqrt or identity link, specified by any of the available ways for GLM links (name, character string, one-element character vector, or object of class \code{link-glm} as returned by \code{\link[stats]{make.link}}). 
}
  \item{trunc}{
  Either \code{0L} for zero-truncated distribution, or \code{-1L} for default untruncated distribution.
  }
}
\details{
The name \code{NB_shape} should be used to set values of shape in control arguments of the fitting functions (e.g., \code{fitme(.,init=list(NB_shape=1))}). 

% For approx as negbin2 with an heteroscedastic Gamma random-effect, specified e.g. as \code{(weights-1|.)} where the weights need to be updated iteratively as function of predicted response, cf file \code{test-negbin1.R} in the \code{/test} directory provides one example. Other mean-variance relationship can be handled in the same way.
}
\value{
A list, formally of class \code{c("LLF", "family")}. See \code{\link{LL-family}} for details about the structure and usage of such objects.
}
\seealso{Examples in \code{\link{LL-family}}. \code{\link{resid.model}} for an example with a residual-dispersion model.}
