\name{negbin}
\alias{negbin}
\alias{negbin2}
\alias{Tnegbin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Family function for negative binomial \dQuote{2} response (including truncated variant).
}
\description{
Returns a GLM \code{\link[stats]{family}} object for negative-binomial model with variance quadratically related to the mean \eqn{\mu}: variance=\eqn{\mu+\mu^2}/shape, with known or unknown underlying Gamma shape parameter. The zero-truncated variant can be specified as \code{negbin2(., trunc = 0L)}. See \code{\link{negbin1}} for the alternative negative-binomial model with variance \dQuote{linearly} related to the mean. 

\code{negbin(.)} is an alias for \code{negbin2(.)} (truncated or not), and \code{Tnegbin(.)} is an alias for \code{negbin2(., trunc = 0L)}.    
}
\usage{
# (the shape parameter is actually not requested unless this is used in a glm() call)
#
negbin2(shape = stop("negbin2's 'shape' must be specified"), link = "log", trunc = -1L, 
        LLgeneric = TRUE)

# For use with glm(), both negbin2's 'shape' and glm's method="llm.fit" are needed.  

# alias with distinct arguments:
Tnegbin(shape = stop("Tnegbin's 'shape' must be specified"), link = "log")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{shape}{
  Shape parameter of the underlying Gamma distribution, given that the \code{negbin} family can be represented as a Poisson-Gamma mixture, where the conditional Poisson mean is \eqn{\mu} times a Gamma random variable with mean 1 and shape \code{shape} (as produced by \code{rgamma(., shape=shape,scale=1/shape)}). 
}
  \item{link}{
  log, sqrt or identity link, specified by any of the available ways for GLM links (name, character string, one-element character vector, or object of class \code{link-glm} as returned by \code{\link[stats]{make.link}}). 
}
  \item{trunc}{
  Either \code{0L} for zero-truncated distribution, or \code{-1L} for default untruncated distribution.
  }
  \item{LLgeneric}{
  For development purposes, not documented. %Boolean: whether to use some generic methods or methods specific for GLM response families. 
  }
}
\details{
\code{shape} is the \eqn{k} parameter of McCullagh and Nelder (1989, p.373) and the \code{theta} parameter of Venables and Ripley (2002, section 7.4). The latent Gamma variable has mean 1 and variance 1/shape. %The \code{negbin} family is sometimes called the NegBin1 model (as the first, historically) in the literature on negative binomial models, and sometimes the NegBin2 model (because its variance is a quadratic function of its mean). 

The name \code{NB_shape} should be used to set values of shape in control arguments of the fitting functions (e.g., \code{fitme(.,init=list(NB_shape=1))}). 
}
\value{
A family object.
}
\references{
McCullagh, P. and Nelder, J.A. (1989) Generalized Linear Models, 2nd edition. London: Chapman & Hall.

Venables, W. N. and Ripley, B. D. (2002) Modern Applied Statistics with S-PLUS. Fourth Edition. Springer.
}
\examples{
## Fitting negative binomial model with estimated scale parameter:
data("scotlip")
fitme(cases~I(prop.ag/10)+offset(log(expec)),family=negbin(), data=scotlip)
negfit <- fitme(I(1+cases)~I(prop.ag/10)+offset(log(expec)),family=Tnegbin(), data=scotlip)
simulate(negfit,nsim=3)
}
\keyword{models}
\keyword{regression}