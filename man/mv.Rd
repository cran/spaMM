\name{mv}
\alias{mv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Virtual factor for multivariate responses
}
\description{
\bold{Motivation}: In a multivariate-response model fitted by \code{fitmv}, one may wish to fit a random-coefficient term appearing in \emph{s} submodels, that is a random effect with realized values for each of these submodels and each group, with values possibly correlated among submodels within groups. Hence one might wish to specify it as a term of the form \code{(<submodel>|group)}, where \code{<submodel>} would represent a factor for the \emph{s} submodels. But the data are not expected to contain a factor for these submodels, so such a syntax would not work without substantial data reshaping. Instead, this effect can be stated as \code{mv(\ldots)} where the \code{\ldots} are the indices of the submodels where the random effect appears. For example if submodels 2 and 3 include this random-coefficient term, the term can be specified as \code{(mv(2,3)|group)}. *No* \dQuote{argument name} should be attached to the indices, as any named element will be treated distinctly (cf. \code{lhs} argument), or ignored.

The \code{mv} syntax can be used further to declare \emph{composite random effects}. As defined in \code{\link{composite-ranef}} and illustrated by a bivariate-response quantitative-genetic model (Examples in \code{\link{Gryphon}}), these effects combine correlations among response variables of different submodels as specified by \code{mv()}, and correlations within each response variable(as specified by a relatedness matrix in the same example).

The \code{mv(\ldots)} expression is treated as a factor for all purposes, meaning or example that \code{(0+mv(2,3)|group)} can also be used, leading (as for any factor) to an alternative parametrization of the same random-coefficient model (see Examples). \bold{The parametrization through \code{(0+mv...)} is generally recommended as its results are easier to interpret}. The random-effect term is treated as a random-coefficient term for all purposes, meaning for example that fixed values can be specified for its parameters using the \code{ranCoefs} syntax (see Examples).  
}
\usage{
# mv(..., xpr)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\ldots}{
  Indices of all the submodels (possibly more than two) where the random effect involving this virtual factor appears. *No* \dQuote{argument name} should be attached to the indices. 
}
  \item{xpr}{
  (experimental) a character string which, when evaluated as an expression in the environment specified by the internal copy of the data, evaluates to a vector (one element per row of the data). \code{as.numeric()} is then applied to the vector, so that the resulting weights are numeric values by which the rows of the design matrix for the random effects will be pre-multiplied. % combined with aliases in some private tests
}
}
\value{
Not a function, hence no return value. In the summary of the fit, levels for the different submodels \emph{s} within each group are labelled \code{.mv}\emph{s}.
}
\seealso{
\code{\link{Gryphon}} for example of composite random effects.

The \code{X2X} argument of \code{\link{fitmv}} for fixed effects shared across sub-models.

}
\examples{
if (spaMM.getOption("example_maxtime")>1.1) {
## data preparation
data("wafers")
me <- fitme(y ~ 1+(1|batch), family=Gamma(log), data=wafers, fixed=list(lambda=0.2))

set.seed(123)
wafers$y1 <- simulate(me, type="marginal")
wafers$y2 <- simulate(me, type="marginal")

## fits
(fitmv1 <- fitmv(
  submodels=list(mod1=list(formula=y1~X1+(mv(1,2)|batch), family=Gamma(log)),
                 mod2=list(formula=y2~X1+(mv(1,2)|batch), family=Gamma(log))), 
  data=wafers))
# alternative '0+' parametrization of the same model:
(fitmv2 <- fitmv(
  submodels=list(mod1=list(formula=y1~X1+(0+mv(1,2)|batch), family=Gamma(log)),
                 mod2=list(formula=y2~X1+(0+mv(1,2)|batch), family=Gamma(log))), 
  data=wafers)) 
# relationship between the *correlated* effects of the two fits
ranef(fitmv2)[[1]][,2]-rowSums(ranef(fitmv1)[[1]]) # ~ 0

# fit with given correlation parameter: 
update(fitmv2, fixed=list(ranCoefs=list("1"=c(NA,-0.5,NA)))) 
}
}
