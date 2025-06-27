\name{pois4mlogit}
\alias{pois4mlogit}
\alias{predict.pois4mlogit}
\alias{reshape2long}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit multinomial logit models.
}
\description{
\code{pois4mlogit} fits a multinomial logit model by wrapping iterative \code{fitmv} fits of a multivariate \code{poisson(log)} surrogate model, according to the following logic. In the (mixed or not) multinomial logit model, the probabilities \eqn{p_{ic}} of the different categories (or \code{types}) \eqn{c=1,...,C} for the \eqn{i}th multinomial draw \eqn{(n_{i1},...,n_{iC})} are of the form 
\deqn{p_{ic}=\frac{e^{\eta_{ic}}}{\sum_{c=1}^C e^{\eta_{ic}}}} 
where each \eqn{\eta_{ic}} is a linear predictor. This model can be fitted using surrogate \code{poisson} models wherein the \code{formula} contains a term \code{offset(.dynoffset)} in addition to the terms specifying each type-specific \eqn{\eta_{.c}}. The procedure will update the \code{.dynoffset} iteratively in the appropriate way to fit the multinomial model.     

Since \code{fitmv} is called, its various specific features can be used, such as random effects correlated across the different response types, or fixed-effect coefficients that may well be shared (by its \code{X2X} argument) across response types. \bold{Note that a shared Intercept term is not identifiable} as implied by the above expression for the \eqn{p_{ic}}s, and non-identifiability can also result when shared predictor values have shared coefficients.

This feature is experimental. The iterative procedure avoids the need for developing specialized software, although it may be slower than such software would be. Not all post-fit procedures may run or return meaningful results from objects returned by \code{pois4mlogit}. There is a specific \code{predict} method for such objects: this method returns by default the predicted frequencies, summing to 1 for each multinomial draw (see Examples).   
}
\usage{
pois4mlogit(submodels, data, to.long=FALSE, init=list(),  control=list(),
  ...,  next_inits=c("ranPars","v_h","fixef"),
  types, n_iter = 1000L, tol=1e-5, fac=1, progress = FALSE)
\method{predict}{pois4mlogit}(object, newdata=NULL, ...)

%reshape2long(mndata, types)
}
\arguments{
  \item{submodels}{
Passed to \code{\link{fitmv}}: see its \code{submodels} argument. This must contain as many \code{poisson(log)} 
submodels as there are response types in the multinomial model. 
}
  \item{data}{
  Data frame; each line contains a multinomial draw and, as usual, the required predictor variables. The counts for the different response types must form different columns of the \code{data} (this is convenient as, e.g., exactly the same data format can be used in binomial fits). The \code{data} are passed unchanged to \code{\link{fitmv}}, unless \code{to.long} is set to TRUE. Any \code{.dynoffset} variable would be internally overwritten, so it generally does not make sense to include such a variable (it would only affect the initial value of the dynamic offset).
}
\item{init}{
list; Passed to \code{fitmv}. Beyond initiating the first \code{fitmv} call, it controls which parameters have explicit initial values in further iterations (though the initial values themselves are then distinct from those of the first iteration).  
}
\item{next_inits}{
Character vector. Controls initiation of a fit from the result of the result of the previous iteration iteration. If \code{"ranPars"} is included, the previous result is used to provide initial values for (\dQuote{outer-optimized}) random-effect parameters. If \code{"v_h"} is included, the previous result is used to provide initial values for the random effects.. If \code{"fixef"} is included, the previous result is used to provide initial values for the fixed effects.
}
\item{to.long}{
Boolean; for optional reformatting of the data (see Details and Examples).
}
\item{control}{
list; passed to \code{\link{fitmv}}.
}
\item{object, newdata}{
\code{pois4mlogit} fit result, and data frame, respectively; passed to \code{\link{predict.HLfit}} (after a local change of its \code{.dynoffset} in the data). 
}
\item{\ldots}{
Further arguments passed to \code{\link{fitmv}} or to \code{\link{predict.HLfit}}. 
}
  \item{types}{
Character vector: type labels, i.e., column names of the counts for the different mutlinomial response categories in \code{data}.
}
  \item{n_iter}{
Integer: maximum number of iteration of iterative algorithm.
}
  \item{tol}{
Numeric: tolerance threshold for determining convergence.
}
  \item{fac}{
Numeric: a control parameter of the iterative algorithm. The default value is better to avoid convergence issues. Higher values may speed up the fit but may also increase the chance of non-convergence, particularly those above 1.5. 
}
  \item{progress}{
Boolean or numeric: whether to print information about number of iterations, and (if \code{progress>1L}) about each iteration. Negative values suppress convergence warnings.
}
}
\details{
At convergence, the dynamic offset should not change over iterations, and the poisson model predictions (offset included) for each multinomial draw should sum to the sample size of the draw. Two convergence criteria, the \code{Ocrit} and the \code{Scrit} respectively, assess these properties. Non-convergence may signal a problem with the fitting procedure (which ca be controlled by its \code{fac} argument), *but* may also signal that the model is not identifiable, in which case it should be modified.

Missing data are handled. This includes the detection and correct handling of cases where (i) a multinomial draw is uninformative because only one type has full information (i.e., response value and all required predictor variables); and (ii) a multinomial draw is informative about some but not all types. In that case, the basic form of the model for a multinomial draw still holds for the relative counts of these types.  

When \code{to.long=TRUE}, the poisson surrogate model is fitted to data specified in a long form where each multinomial draw of size \eqn{s_i} is described as \eqn{s_i} multinomial draws of size 1 (so that each Poisson submodel is fitted to a response vector of 0s and 1s). In this long format, counts for the different response types still form different columns (here containing only 0 or 1) of the \code{data}. This format is not needed (and not memory-efficient) but it may be useful to compare more easily the likelihood of the surrogate model to the multinomial one (see Examples), and it was used in such a way by Chen & Kuo (2001). The Examples use the utility function \code{reshape2long} to convert the original data frame to the long format. 
}
\value{
A \code{list} also inheriting from class \code{HLfit} (as the the return value of the \code{fitmv} call from which the present value is derived), and from class \code{pois4mlogit}.
}
\references{
The use of an iteratively updated offset is inspired from the arguments in Chen & Kuo (2001), who described related approaches to fit multinomial logit models using poisson(log) surrogate ones.

Chen, Z. and Kuo, L. (2001) A note on the estimation of the multinomial logit model with random
effects. The American Statistician 55, 89-95. \url{https://www.jstor.org/stable/2685993}
}
\examples{
#### Fitting a binomial(logit) model by a bivariate poisson(log) surrogate:
## Toy data: Let us say we observe a color polymorphism of irises in 10 populations...
set.seed(123)
ssize <- 10L
shape <- 0.35
toydata <- data.frame(
  yellow=rbinom(ssize, 16, prob=rbeta(ssize,shape,shape)),
  purple=rbinom(ssize, 16, prob=rbeta(ssize,shape,shape)), # (purple ignored below)
  blue=rbinom(ssize, 16, prob=rbeta(ssize,shape,shape)),
  phenotype=rnorm(ssize)
)
## Standard binomial fit (purple flowers are ignored here)
(byB <- fitme(cbind(yellow,blue) ~ phenotype, family = binomial(), 
      data=toydata))

## Surrogate fit
(byP2 <- pois4mlogit(submodels = list(
  list(yellow ~ offset(.dynoffset) + phenotype, family = poisson()),
  list(blue ~ offset(.dynoffset) + phenotype, family = poisson())),
  data = toydata, types=c("yellow","blue")))
  
#### Recover summaries of the binomial fit from the surrogate fit:
  
## Coefficients of binomial model recovered as:
  fixef(byP2)[1:2]-fixef(byP2)[3:4]
  
## SEs of coefficients of binomial model recovered as:  
  {
    P2B <- rbind(c(1,0,-1,0),c(0,1,0,-1))
    P2B \%*\% vcov(byP2) \%*\% t(P2B)
  } 
  # practically equivalent to 
  vcov(byB)
  
## Probabilities of each type:
  matrix(predict(byP2),ncol=2,
         dimnames=list(NULL, c("yellow","blue")))           
  
## logLiks
# The logLiks differ between the binomial and poisson surrogate fits 
# because the combinatorial coefficients differ. The relationship
# between these logLiks is simple when the long form of the data is used:
  str(long2 <- reshape2long(toydata, c("yellow","blue")))
  
# Fits on long data:
  (byBlong  <- fitme(cbind(yellow,blue) ~ phenotype, family = binomial(), 
                     data=long2))
  (byPlong <- pois4mlogit(submodels = list(
    list(yellow ~ offset(.dynoffset) + phenotype, family = poisson()),
    list(blue ~ offset(.dynoffset) + phenotype, family = poisson())),
    data = toydata, to.long=TRUE, types=c("yellow","blue")))
    
# The logLik of the surrogate model is
  sum(byPlong$eta*byPlong$y) - sum(exp(byPlong$eta)) # = logLik(byPlong)
# while le logLik of the binomial one can be obtained as 
  sum(byPlong$eta*byPlong$y) # = logLik(byBlong)
# the difference is 156  __= total multinomial sample size__
    
# This illustrates that the logLiks of different surrogate models differ 
# only by a constant independent of the fitted model. 
# Likelihood ratios between surrogate models are then 
# identical to likelihood ratios between corresponding binomial models,
# provided a single data format is used throughout the comparisons.

}

