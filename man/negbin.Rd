\name{negbin}
\alias{negbin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Family function for negative binomial GLMs and mixed models.
}
\description{
Specifies the information required to fit a negative binomial generalized linear model, 
with known or unknown underlying Gamma shape parameter.
}
\usage{
negbin(shape = stop("negbin's 'shape' must be specified"), link = "log")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{shape}{
  Shape parameter of the underlying Gamma distribution, given that the \code{negbin} family can be represented as a Poisson-Gamma mixture, where the conditional Poisson mean is \eqn{\mu} times a Gamma random variable with mean 1 and shape \code{shape} (as produced by \code{rgamma(., shape=shape,scale=1/shape)}). 
}
  \item{link}{
  log, sqrt or identity link, specified by the sevaral available ways for GLM links (name, character string, one-element character vector, or object of class \code{link-glm} as returned by \code{\link[stats]{make.link}}). 
}
}
\details{
\code{shape} is the \eqn{k} parameter of McCullagh and Nelder (1989, p.373) and the \code{theta} parameter of Venables and Ripley (2002, section 7.4). The \code{negbin} family is sometimes called the NegBin1 model in the literature on negative binomial models.
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
data(scotlip)
fitme(cases~I(prop.ag/10)+offset(log(expec)),family=negbin(), data=scotlip)
}
\keyword{models}
\keyword{regression}