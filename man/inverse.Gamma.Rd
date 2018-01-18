\name{inverse.Gamma}
\alias{inverse.Gamma}
\alias{Gamma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Distribution families for Gamma and inverse Gamma-distributed random effects
}
\description{
For dispersion parameter \eqn{\lambda}, \code{Gamma} means that random effects are distributed as \eqn{u ~}Gamma(\code{shape=}1/\eqn{\lambda},\bold{scale=}\eqn{\lambda}), so \eqn{u} has mean 1 and variance \eqn{\lambda}. Both the log (\eqn{v=log(u)}) and identity (\eqn{v=u}) links are possible, though in the latter case the variance of \eqn{u} is constrained below 1 (otherwise Laplace approximations fail).

The two-parameter inverse Gamma distribution is the distribution of the reciprocal of a variable distributed according to the Gamma distribution Gamma with the same shape and scale parameters. \code{inverse.Gamma} implements the one-parameter inverse Gamma family with shape=1+1/\eqn{\lambda} and \bold{rate}=1/\eqn{\lambda}) (rate=1/scale). It is used to model the distribution of random effects. Its mean=1; and its variance =\eqn{\lambda/(1-\lambda))} if \eqn{\lambda<1}, otherwise infinite. The default link is \code{"-1/mu"}, in which case \code{v=-1/u} is \dQuote{-Gamma}-distributed with the same shape and rate, hence with mean \eqn{-(\lambda+1)} and variance \eqn{\lambda(\lambda+1)}, which is a different one-parameter Gamma family than the above-described \code{Gamma}. The other possible link is \code{v=log(u)} in which case\cr 
\eqn{v ~ -log(X~}Gamma\eqn{(1+1/\lambda,1/\lambda))}, with mean \eqn{-(log(1/\lambda)+}digamma\eqn{(1+1/\lambda))} and variance trigamma(\eqn{1+1/\lambda}).
}
\usage{
inverse.Gamma(link = "-1/mu")
# Gamma(link = "inverse") using stats::Gamma 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{link}{For \code{Gamma}, allowed links are \code{log} and \code{identity} (the default link from \code{\link[stats]{Gamma}}, \code{"inverse"}, cannot be used for the random effect specification). For \code{inverse.Gamma}, allowed links are \code{"-1/mu"} (default) and \code{log}.}
}
\examples{
# see help("HLfit") for fits using the inverse.Gamma distribution.
}
