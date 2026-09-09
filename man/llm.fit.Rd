\name{llm.fit}
\alias{LL-family}
\alias{llm.fit}
\title{
Link-linear regression models (LLMs)
}
\description{
Some \dQuote{family} objects in \pkg{spaMM} describe models with non-GLM distribution families, such as the \code{\link{negbin1}} or \code{\link{beta_resp}} families already widely considered in previous works and other packages. These models are characterized by a linear predictor, a link function, and a distribution for residual variation that does not belong to the exponential family from which GLMs are defined. 

These family objects are conceived for use with \pkg{spaMM}'s fitting functions. They cannot generally be used as argument to the \code{glm} function, except when this function is highjacked by use of the \code{method="llm.fit"} argument, where llm stands for Link-Linear (as in \dQuote{log-linear}, say) regression Model. 

Mixed-effect models fitted by such methods cannot use approximations based on the expected information matrix, in contrast to GLM distribution families. 
}
\usage{
# glm(..., method="llm.fit")
## See also 'beta_resp', 'negbin1', 'betabin', and possibly later additions.
}
\details{
These family objects are lists, formally of class \code{c("LLF", "family")}. Compared to a \code{\link[stats]{family}} object, they have additional elements, not documented here. 

As \code{stats::} GLM family objects do, they provide deviance residuals through the \code{dev.resids} member function. There are various definitions of deviance residuals for non-GLM families in the literature. Here they are defined from the unit deviances \dQuote{2*(saturated_logLik - logLik)\eqn{_i}} for each observation \eqn{i} (provided by the \code{dev.resids} function defined by the family), where the likelihood for the saturated model is the likelihood maximized wrt to the mean parameter \eqn{\mu_i} for each observation \eqn{y_i} independently. The maximizing  \eqn{\mu_i} is not equal to the observation \eqn{y_i}, in contrast to the standard result for GLMs.   

\code{llm.fit} provides the deviance (the sum of the above unit deviances), but not the null deviance provided by \code{stats::glm.fit}. The latter is the deviance of an intercept-only (or offset only) model. There, the intercept estimate is the mean response value for GLMs, so that the null deviance is immediately obtained; but this is not generally so for other families. 
}
\examples{
data(scotlip)

### negbin1 response:

# Fixed-effect model
#
(var_shape <- fitme(cases~I(prop.ag/10)+offset(log(expec)),family=negbin1(), 
                    data=scotlip))

# Highjacking glm(): the family parameter must be given 
#
fitted_shape <- residVar(var_shape,which="fam_parm")
glm(cases~I(prop.ag/10)+offset(log(expec)),family=negbin1(shape=fitted_shape), 
    method="llm.fit", data=scotlip)

### Similar exercice with Beta family: 

set.seed(123)
beta_dat <- data.frame(y=runif(100),grp=sample(2,100,replace = TRUE))
  
# Fixed-effect model
(var_prec <- fitme(y ~1, family=beta_resp(), data= beta_dat))
  
# Highjacking glm(): 
fitted_prec <- residVar(var_prec,which="fam_parm")
glm(y ~1, family=beta_resp(prec=fitted_prec), data= beta_dat, method="llm.fit") 

}
