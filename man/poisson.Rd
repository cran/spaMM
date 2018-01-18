\name{Poisson}
\alias{Poisson}
\alias{Tpoisson}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Family function for GLMs and mixed models with Poisson and zero-truncated Poisson response.
}
\description{
\code{Poisson} (with a capital P) is a \code{\link[stats]{family}} that specifies the information required to fit a Poisson generalized linear model. Differs from the base version \code{stats::poisson} only in that it handles the zero-truncated variant, which can be specified either as \code{Tpoisson(<link>)} or as \code{Poisson(<link>, trunc = 0L)}.
}
\usage{
Poisson(link = "log", trunc = -1L)
Tpoisson(link="log")
}
\arguments{
  \item{link}{
  log, sqrt or identity link, specified by the sevaral available ways for GLM links (name, character string, one-element character vector, or object of class \code{link-glm} as returned by \code{\link[stats]{make.link}}). 
}
  \item{trunc}{
  Either \code{0L} for zero-truncated distribution, or \code{-1L} for default untruncated dsitribution.
  }
}
\value{
A family object.
}
\references{
McCullagh, P. and Nelder, J.A. (1989) Generalized Linear Models, 2nd edition. London: Chapman & Hall.
}
\examples{
data("scotlip")
logLik(glm(I(1+cases)~1,family=Tpoisson(),data=scotlip))
logLik(fitme(I(1+cases)~1+(1|id),family=Tpoisson(),fixed=list(lambda=1e-8),data=scotlip))
}
\keyword{models}
\keyword{regression}