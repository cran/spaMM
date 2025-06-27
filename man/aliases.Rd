\name{aliases}
\alias{aliases}

\title{Variable aliases for multivariate-response fits}
\description{
Consider the \code{toy_iris} dataset defined in the Examples, noting in particular
\preformatted{
u <- rnorm(10,sd=1)
id_y <- gl(10,5)
id_b <- rep(seq(10),5)
}
an the subsequent use of \code{u[id_y]} and \code{u[id_b]} in the expectations of Poisson draws. Suppose that we try to a multivariate-response model with two Poisson responses to it (yes, it is an artificial example, as a multivariate fit may not be necessary here, but if you know what \code{\link{pois4mlogit}} does, then a less artificial example can be defined from a similar dataset with binomial samples instead of Poisson ones). A fit by 
\preformatted{
fitmv(submodels = list(
    list(yellow ~ 1+(1|id_y), family = poisson()),
    list(blue ~ 1+(1|id_b), family = poisson())),
    data = toy_iris)
}
does not match the data-generating algorithm, because the fit of the \code{(1|id_y)} and \code{(1|id_b)} random effects does not take into account that the latent values of these random effects are sampled from a single vector \code{u} of 10 latent values. A matching fit should fit a single random-effect variance and produce a single predicted vector for \code{u}, rather than two distinct vectors of predicted values. 

Conversely, a fit such as  
\preformatted{
fitmv(submodels = list(
    list(yellow ~ 1+(1|id_y), family = poisson()),
    list(blue ~ 1+(1|id_y), family = poisson())),
    data = toy_iris)
}
would fit a single random-effect variance and produce a single predicted vector, but the factor indices are incorrect in the second submodel.

We need a syntax such that a single variance and a single vector are fitted (as when the random-effect terms are identical in the two formulas), but where the factor indices are effectively \code{id_y} and \code{id_b} in the two submodels. The syntax of the \code{proper} fit in the Examples achieves this effect. The random effect is specified as 
\code{(1|id)} where \code{id} is *not* an actual factor in the \code{data} but is an alias for the variable \code{id_y} in the first formula and \code{id_b} in the second one. This interpretation of the \code{(1|id)} term is specified by the argument
\code{aliases=list(id=c("id_y","id_b"))}.
}
\examples{
# Toy data
set.seed(123)
ssize <- 50L
u <- rnorm(10,sd=1)
id_y <- gl(10,5)
id_b <- rep(seq(10),5)
yellow <- rpois(ssize, lambda=exp(1+u[id_y]))
blue <- rpois(ssize, lambda=exp(1+u[id_b]))
toy_iris <- data.frame(
  id_y=id_y, id_b=id_b, yellow=yellow, blue=blue
)

# Fit
proper <- fitmv(submodels = list(
    list(yellow ~ 1+(1|id), family = poisson()),
    list(blue ~ 1+(1|id), family = poisson())),
    data = toy_iris, aliases=list(id=c("id_y","id_b")))
    
ranef(proper) # single vector of 10 values

# Rows 1-50 and 51-100 of Z matrix show 
# the distinct design of the two submodels:
get_matrix(proper,"ZA")
}


