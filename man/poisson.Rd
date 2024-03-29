\name{Poisson}
\alias{Poisson}
\alias{Tpoisson}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Family function for GLMs and mixed models with Poisson and zero-truncated Poisson response.
}
\description{
\code{Poisson} (with a capital P) is a \code{\link[stats]{family}} that specifies the information required to fit a Poisson generalized linear model. Differs from the base version \code{stats::poisson} only in that it handles the zero-truncated variant, which can be specified either as \code{Tpoisson(<link>)} or as \code{Poisson(<link>, trunc = 0L)}. The truncated poisson with mean \eqn{\mu_T} is defined from the un-truncated poisson with mean \eqn{\mu_U}, by restricting its response strictly positive value. \eqn{\mu_T=\mu_U/(1-p0)}, where \eqn{p0:=\exp(-\mu_U)} is the probability that the response is 0. 
}
\usage{
Poisson(link = "log", trunc = -1L, LLgeneric=TRUE)
Tpoisson(link="log")
# <Poisson object>$linkfun(mu, mu_truncated = FALSE)
# <Poisson object>$linkinv(eta, mu_truncated = FALSE)
}
\arguments{
  \item{link}{
  log, sqrt or identity link, specified by any of the available ways for GLM links (name, character string, one-element character vector, or object of class \code{link-glm} as returned by \code{\link[stats]{make.link}}). 
}
  \item{trunc}{
  Either \code{0L} for zero-truncated distribution, or \code{-1L} for default untruncated distribution.
  }
  \item{eta,mu}{Numeric (scalar or array). The linear predictor; and the expectation of response, truncated or not depending on \code{mu_truncated} argument. }
  \item{mu_truncated}{Boolean. For \code{linkinv}, whether to return the expectation of truncated (\eqn{\mu_T}) or un-truncated (\eqn{\mu_U}) response. For \code{linkfun}, whether the \code{mu} argument is \eqn{\mu_T}, or is \eqn{\mu_U}  but has \eqn{\mu_T}  as attribute (\eqn{\mu_U} without the attribute is not sufficient). }
    \item{LLgeneric}{
  For development purposes, not documented. %Boolean: whether to use some generic methods or methods specific for GLM response families. 
  }

}
\value{
A \code{family} object suitable for use with \code{glm}, as \code{stats::} family objects.
}
\details{
Molas & Lesaffre (2010) developed expressions for deviance residuals for the truncated Poisson distribution, which were the ones implemented in \pkg{spaMM} until version 3.12.0. Later versions implement the (non-equivalent) definition as \dQuote{2*(saturated_logLik - logLik)}.   

\code{predict}, when applied on an object with a truncated-response family, by default returns \eqn{\mu_T}. The simplest way to predict \eqn{\mu_U} is to get the linear predictor value by \code{predict(.,type="link")}, and deduce \eqn{\mu_U} using \code{linkinv(.)} (with default argument \code{mu_truncated=FALSE}), since getting \eqn{\mu_U} from \eqn{\mu_T} is comparatively less straightforward. The \code{mu.eta} member function is that of the base \code{poisson} family, hence its \code{mu} argument represents \eqn{\mu_U}. 

\code{simulate}, when applied on an object with a truncated-response family, simulates the truncated family. There is currently no clean way to override this (trying to pass\code{type="link"} to \code{predict} will not have the intended effect). % wrong predict type for mu + no effect on r_resid_var conditional simulation.   
}
\references{
McCullagh, P. and Nelder, J.A. (1989) Generalized Linear Models, 2nd edition. London: Chapman & Hall.

Molas M. and Lesaffre E. (2010). Hurdle models for multilevel zero-inflated data via h-likelihood. Statistics in Medicine 29: 3294-3310.
}
\examples{
data("scotlip")
logLik(glm(I(1+cases)~1,family=Tpoisson(),data=scotlip))
logLik(fitme(I(1+cases)~1+(1|id),family=Tpoisson(),fixed=list(lambda=1e-8),data=scotlip))
}
\keyword{models}
\keyword{regression}