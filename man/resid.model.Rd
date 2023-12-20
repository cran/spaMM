\name{resid.model}
\alias{resid.model}

\title{Structured dispersion models}
\description{
The \code{resid.model} argument of fitting functions can be used to specify a model for a residual-dispersion parameter of various response families, that is, either\cr
(1) the \eqn{\phi} parameter of the gaussian and Gamma GLM families;\cr
(2a) the dispersion parameter of some other GLM families, such as the shape parameter of the negbin1 and negbin2 families; or\cr
(2b) the dispersion parameter of some other (non-GLM) response families, such as the precision parameter of the beta response family.

This documentation is more specifically for case (2). Case (1) is more specifically documented as \code{\link{phi-resid.model}}.

In case (2) the model for the dispersion parameter is constrained as a fixed-effect model, of the form\cr
dispersion parameter =\code{ exp(}\bold{X}\eqn{\beta}\code{+offset}),\cr 
and specified using the standard \code{formula} syntax. Random effects cannot be included, in contrast to dispersion models for case (1).
}
\usage{
# 'resid.model' argument of fitme() and fitmv()
}
\arguments{
The \code{resid.model} for case (2) is simply a formula (without left-hand side) for the logarithm of the dispersion parameter. Fixed \eqn{\beta} values can be specified through the \code{rdisPars} element of the \code{fixed} argument in the \code{fitme} call (or through the \code{fixed} argument of each submodel of a \code{fitmv} call). Likewise, initial values can be specified through the \code{init} argument.  
}
\details{
In case (2) a fixed \dQuote{heteroscedastic} model can also be specified directly through the family specification, e.g., family=negbin1(shape=<vector>) where the vector has the length of the response vector, but this may not be suitable if the model is to be used for prediction purposes (where the residual-dispersion model should be specified in such a way that one can \dQuote{predict} new dispersion values from it).

The design matrix for the specified model is internally rescaled to avoid numerical problems. That means that there is no need to rescale the predictor variable, even if it tends to take large (cf \sQuote{population} variable in the Examples) of small values (this is also true for fixed-effect predictors of the mean-response model).
}
\value{
The fit results for the residual model are accessible through the summary and various extractors. In particular, the \code{get_fittedPars} extractor will by default include in its return value the \code{rdisPars} element, which is here the vector of fitted \eqn{\beta} coefficients.  \code{residVar(., which="fam_parm")} will return the vector of fitted values of the dispersion parameter. 
}
\examples{
data("scotlip")
if (spaMM.getOption("example_maxtime")>3) {
(toyfit <- fitme(cases~1+(1|id),family=negbin1(), data=scotlip, resid.model = ~ population))

# => This toy example is a bit challenging to fit because the data set is small and 
# individual-level variation is here described both by a random effect 
# and by a two-parameter negbin1 residual variation. The fit might often stop 
# at a local maximum of the logLik in such cases (although there is no evidence
# that this is presently the case).
}
}