\name{beta_resp}
\alias{beta_resp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Beta-response family object
}
\description{
Returns a \code{family} object for beta-response models.
The model described by such a family is characterized by a linear predictor, a link function, and the beta density for the residual variation.

The precision parameter \code{prec} of this family is a positive value such that the variance of the response given its mean \eqn{\mu} is \eqn{\mu(1-\mu)/(1+}\code{prec}). \code{prec} is thus the precision parameter \eqn{\phi} of Ferrari & Cribari-Neto (2004) and of the \pkg{betareg} package (Cribari-Neto & Zeileis 2010).  

Prior weights are meaningful for this family and handled as a factor of the precision parameter (as for GLM families) hence here not as a divisor of the variance (in contrast to GLM families): the variance of the response become \eqn{\mu(1-\mu)/(1+}\code{prec*<prior weigths>}).
}
\usage{
beta_resp(prec = stop("beta_resp's 'prec' must be specified"), link = "logit")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prec}{
  Scalar (or left unspecified): precision parameter of the beta distribution. 
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
  beta_dat <- data.frame(y=runif(100),grp=sample(2,100,replace = TRUE))
  
  fitme(y ~1+(1|grp), family=beta_resp(), data= beta_dat)
  ## same logL, halved 'prec' when prior weights=2 are used: 
  # fitme(y ~1+(1|grp), family=beta_resp(), data= beta_dat, prior.weights=rep(2,100))
}
