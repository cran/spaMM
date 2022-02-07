\name{residVar}
\alias{residVar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Residual variance extractor
}
\description{
Extracts from a fit object the residual variance or, depending on the \code{which} argument, a family dispersion parameter, or a vector of values of the dispersion parameter phi (which is not the residual variance itself for gamma-response models), or further information about the residual variance model. 
}
\usage{
residVar(object, which = "var", submodel = NULL, newdata = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of class \code{HLfit}, as returned by the fitting functions in \code{spaMM}.}
  \item{which}{Character: \code{"var"} for the fitted residual variances, \code{"phi"} for the fitted phi values, \code{"fam_parm"} for the dispersion parameter of \code{COMPoisson} and \code{negbin} families, \code{"fit"} for the fitted residual model (a GLM or a mixed model for residual variances, if not a simpler object),  and \code{"family"} or \code{"formula"} for such properties of the residual model.}
  \item{submodel}{integer: the index of a submodel, if \code{object} is a multivariate-response model fitted by \code{fitmv}. This argument is mandatory for all \code{which} values except \code{"var"} and \code{"phi"}.}
    \item{newdata}{\bold{Either} NULL, a matrix or data frame, or a numeric vector. See \code{\link{predict.HLfit}} for details.}
}
\value{
Default \code{which="var"} (or \code{"phi"}) always return a vector of residual variances (or, alternatively, phi values) of length the number of responses of the fit (or only the number of responses for a submodel, if relevant). \code{which="fit"} returns an object of class \code{HLfit}, \code{glm}, or a single scalar depending on the residual dispersion model. Other \code{which} values return an object of class \code{family} or \code{formula} as expected. 
}
\seealso{
\code{\link{get_residVar}} is a alternative extractor of residual variances with different features inherited from \code{get_predVar}. In particular, it is more suited for computing the residual variances of new realizations of a fitted model, not accounting for prior weights used in fitting the model (basic examples of using the \code{IsoriX} package provide a context where this is the appropriate design decision). By contrast, \code{residVar} aims to account for prior weights.
}
\examples{
# data preparation: simulated trivial life-history data
set.seed(123)
nind <- 20L
u <- rnorm(nind)
lfh <- data.frame(
  id=seq_len(nind), id2=seq_len(nind), 
  feco= rpois(nind, lambda = exp(1+u)), 
  growth=rgamma(nind,shape=1/0.2, scale=0.2*exp(1+u)) # mean=exp(1+u), var= 0.2*mean^2
)
# multivariate-response fit                  
fitlfh <- fitmv(submodels=list(list(feco ~ 1+(1|id), family=poisson()),
                               list(growth ~ 1+(1|id), family=Gamma(log))),
                data=lfh)
#
residVar(fitlfh)
residVar(fitlfh, which="phi") # shows fixed phi=1 for Poisson responses
residVar(fitlfh, submodel=2)
residVar(fitlfh, which="family", submodel=2)
residVar(fitlfh, which="formula", submodel=2)
residVar(fitlfh, which="fit", submodel=2) # Fit here characterized by a single scalar
}
