\name{DoF}
\alias{DoF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Degrees of freedom extractor
}
\description{
This extracts the number of degrees of freedom for a model, in the usual sense for likelihood-ratio tests: a count of number of fitted parameters, distinguishing different classes of parameters (see Value).
}
\usage{
DoF(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
A fitted-model object, of class \code{"HLfit"}.
}
}
\details{
The output distinguishes counts of random-effect vs residual-dispersion parameters, following the conceptual distinction between effects that induce correlations between different levels of the resonse vs. observation-level effects. However, a residual-dispersion component can be declared as a random effect, so that the counts for logically equivalent models may differ according to the way a model was declared. For example if residual dispersion for an LLM is declared as an observation-level random effect while \code{phi} is fixed, the \code{p_lambda} component will include 1 df for what would otherwise be accounted by the \code{p_rdisp} component. A more involved case where the same contrast happens is when a negative-binomial model (with a residual-dispersion \code{shape} parameter) is declared as a Poisson-gamma mixture model (with a varaince parameter for the Gamma-distributed individual-level random effect).    
}
\value{
A vector with possible elements \code{p_fixef}, \code{p_lambda}, \code{p_corrPars} and \code{p_rdisp} for, respectively, the number of fixed-effect coefficients of the main-response model, the number of random-effect variance parameters, the number of random-effect correlation parameters, and the number of residual dispersion parameters (the latter being itself, for a mixed-effect residual-dispersion model, the sum of such components).
}
\seealso{
\code{\link{df.residual.HLfit}}; \code{\link{get_any_IC}} for extracting effective degrees of freedom considered in the model-selection literature; \code{\link{as_LMLT}} for access to the effective degrees of freedom considered in Satterthwaite's test and its extentions.
}

