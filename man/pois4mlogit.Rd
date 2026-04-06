\name{pois4mlogit}
\alias{pois4mlogit}
\alias{predict.pois4mlogit}
\alias{reshape2long}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit multinomial logit models and multivariate-response models including them.
}
\description{
\code{pois4mlogit} is a procedure that fits a multinomial logit model by calling \code{fitmv} to fit multivariate \code{poisson(log)} surrogate models, according to the following logic. In the (mixed or not) multinomial logit model, the probabilities \eqn{p_{ic}} of the different categories (or \code{types}) \eqn{c=1,...,C} for the \eqn{i}th multinomial draw \eqn{(n_{i1},...,n_{iC})} are of the form 
\deqn{p_{ic}=\frac{e^{\eta_{ic}}}{\sum_{c=1}^C e^{\eta_{ic}}}} 
where each \eqn{\eta_{ic}} is a linear predictor. By contrast, the denominator makes \eqn{\log(p_{ic})} non-linear. This type of model can can be fitted as a Poisson non-linear mixed-effect model (e;g., Chen & Kuo, 2001). \pkg{spaMM} does not have general procedures for fitting non-linear mixed-effect models, but its current procedures for GLMMs have been hacked to fit this specific model, by combining (1) the iterative adjustment of GLMMs where the \eqn{\log(\sum_{c=1}^C e^{\eta_{ic}})} values are represented by an offset term, modified over iterations; and (2) ad-hoc modifications of the gradient vector and information matrix of the GLMMs. The correctness of the procedure can be verified for binomial data by comparison to a standard binomial GLMM fit. 

In a \code{pois4mlogit} call, the \code{formula} for each type of the multinomial response must contain a term \code{offset(.dynoffset)} that declares the above-defined offset term, in addition to the terms specifying each type-specific \eqn{\eta_{.c}}. The values of this offset for all \eqn{i} are iteratively updated between successive calls to surrogate Poisson multivariate-response GLMMs, each fitted using \code{fitmv}. 

This procedure is \bold{experimental}. Not all post-fit procedures may run or return meaningful results from objects returned by \code{pois4mlogit}. There is a specific \code{predict} method for such objects: this method returns by default the predicted frequencies, summing to 1 for each multinomial draw (see Examples). \code{\link{plot_effects}} and \code{\link{LRT}} can also handle \code{pois4mlogit} fits.

In mixed-effect models, the procedure iteratively updating the ad-hoc offset has to be called for each combination of random-effect parameters considered by the overall fitting procedure. This may be much slower than more specialized software would be, an issue that may be negligible in some applications and bothering in others. This issue is left for possible later improvements.   

Since \code{fitmv} is called, its various specific features can be used, such as random effects correlated across the different response types (\code{\link{mv}} syntax), or the additional relationships between random effects that can be specified by the \code{\link{aliases}} argument, or shared fixed-effect coefficients across response types (\code{\link{X2X}} argument). Further, additional submodels not belonging to the multinomial model can also be included. However, \code{pois4mlogit} cannot fit a model including several multinomial models.

Users should quickly learn to avoid \bold{non-identifiable models}, in particular the following two simple cases: (i) models where all submodels have their own intercept; and (ii) models where a fixed-effect term (which may be simply an intercept) has the same value in the different submodels, as occurs if its fitted coefficient is shared among submodels, and the values of its regressor are also identical among submodels.
%as implied by the above expression for the \eqn{p_{ic}}s.

}
\usage{
pois4mlogit(submodels, data, to.long=FALSE, init=list(),  
            control=list(wdfac=2, p4m="", grad=FALSE), ...,  
            next_inits=c("ranPars","v_h","fixef"),
            types, n_iter = 1000L, tol=c(1e-3,1e-5), 
            initfn=get_inits_from_fit, progress = FALSE)
\method{predict}{pois4mlogit}(object, newdata=NULL, ...)

reshape2long(data, types) 
}
\arguments{
  \item{submodels}{
Passed to \code{\link{fitmv}}: see its \code{submodels} argument. This must contain as many \code{poisson(log)} 
submodels as there are response types in the multinomial model. 
}
  \item{data}{
  Data frame; each line contains a multinomial draw and, as usual, the required predictor variables. The counts for the different response types must form different columns of the \code{data} (this is convenient as, e.g., exactly the same data format can be used in binomial fits). The \code{data} are passed unchanged to \code{\link{fitmv}}, unless \code{to.long} is set to TRUE. An initial \code{.dynoffset} can be provided in the \code{data}. It will be internally overwritten after the first iteration.
}
\item{init}{
list; Passed to \code{fitmv}. Beyond initiating the first \code{fitmv} call, it controls which parameters have explicit initial values in further iterations (though the initial values themselves are then distinct from those of the first iteration).  
}
\item{next_inits}{
Character vector. Controls initiation of the next \code{fitmv} call from the result of the previous one. If \code{"ranPars"} is included, the function specified by the \code{initfn} argument is applied to this previous fit to provide initial values for (\dQuote{outer-optimized}) random-effect parameters in the next \code{fitmv} call. See Details for possible alternative to the \emph{next-init default}.
If \code{"v_h"} is included, the previous \code{fitmv} result is used to provide initial values for the random effects. If \code{"fixef"} is included, the previous \code{fitmv} result is used to provide initial values for the fixed effects.
}
\item{to.long}{
Boolean; for optional reformatting of the data (see Details and Examples).
}
\item{control}{
list; passed to \code{\link{fitmv}}, but may also contain extra elements intrinsic to \code{pois4mlogit}, as shown in its default value. These elements are currently for development purposes. 
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
Numeric: tolerance thresholds for determining convergence. For development purposes, not documented.
}
%  \item{fac}{
%Numeric: a control parameter of the iterative algorithm. The default value is better to avoid convergence issues. Higher values may speed up the fit but may also increase the chance of non-convergence, particularly those above 1.5. 
%}
  \item{initfn}{
function: a function that returns initial values. Currently for development purposes. % ____F I X M E____
}
  \item{progress}{
Boolean or numeric: whether to print information about number of iterations, and (if \code{progress>1L}) about each iteration. Negative values suppress convergence warnings.
}
}
\details{
At convergence, the dynamic offset should not change over iterations, and the poisson model predictions (offset included) for each multinomial draw should sum to the sample size of the draw. Two convergence criteria, the \code{Ocrit} and the \code{Scrit} respectively, assess these properties. Non-convergence may signal a problem with the fitting procedure (which ca be controlled by its \code{fac} argument), *but* may also signal that the model is not identifiable, in which case it should be modified.

Missing data are handled. This includes the detection and correct handling of cases where (i) a multinomial draw is uninformative because only one type has full information (i.e., response value and all required predictor variables); and (ii) a multinomial draw is informative about some but not all types. In that case, the basic form of the model for a multinomial draw still holds for the relative counts of these types.  

When \code{to.long=TRUE}, the poisson surrogate model is fitted to data specified in a long form where each multinomial draw of size \eqn{s_i} is described as \eqn{s_i} multinomial draws of size 1 (so that each Poisson submodel is fitted to a response vector of 0s and 1s). In this long format, counts for the different response types still form different columns (here containing only 0 or 1) of the \code{data}. This format is not needed (and not memory-efficient) but it may be useful to compare more easily the likelihood of the surrogate model to the multinomial one (see Examples), and it was used in such a way by Chen & Kuo (2001). The Examples use the utility function \code{reshape2long} to convert the original data frame to the long format. 

Some control of initial values for (\dQuote{outer-optimized}) random-effect parameters over iterations is possible. However, the following non-default options are presumably very inefficient and should be considered only when the default option fails to provide appropriate results. 
The \emph{next-init default}, defined in the documentation of the \code{next_inits} argument, can be modified as follows: if \code{"ranPars"} is absent but \code{"init"} is present in the \code{next_inits} vector, this default is modified by the user-provided \code{init} value. If both \code{"ranPars"} and \code{"init"} are absent from the \code{next_inits} vector, user-provided \code{init} values are not used but their names define the elements from the next-init default that are retained in the initial values of the next \code{fitmv} call (thus, other elements are dropped and their initial values will \bold{not} be given by results of the previous \code{fitmv}).
}
\value{
A \code{list} also inheriting from class \code{HLfit} (as the the return value of the \code{fitmv} call from which the present value is derived), and from class \code{pois4mlogit}.
}
\references{
Chen, Z. and Kuo, L. (2001) A note on the estimation of the multinomial logit model with random
effects. The American Statistician 55, 89-95. \url{https://www.jstor.org/stable/2685993}
}
\seealso{
Specific documentations for \code{fitmv} arguments handled by \code{pois4mlogit}, describing relationships between model terms from different sub-models:
\code{\link{mv}} and \code{\link{aliases}} syntaxes for random effects, and \code{\link{X2X}} for fixed effects. 
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
  list(blue ~ offset(.dynoffset) + 0, family = poisson())),
  data = toydata, types=c("yellow","blue")))
  
# => Note the different standard errors: the poisson surrogate does not
# provide the "correct" standard errors of the binomial fit.
# The 'over-parametrized model' below gives standard errors
# closer to the ones from the standard binomial fit.

# Add a trivial gaussian submodel just to show that this can be done:   
(just4fun <- pois4mlogit(submodels = list(
  list(yellow ~ offset(.dynoffset) + phenotype, family = poisson()),
  list(blue ~ offset(.dynoffset) + 0, family = poisson()),
  list(phenotype ~1)),
  data = toydata, types=c("yellow","blue")))

#### Over-parametrized model:

(byP2over <- pois4mlogit(submodels = list(
  list(yellow ~ offset(.dynoffset) + phenotype, family = poisson()),
  list(blue ~ offset(.dynoffset) + phenotype, family = poisson())),
  data = toydata, types=c("yellow","blue")))
  
## Coefficients of binomial model recovered as:
  fixef(byP2over)[1:2]-fixef(byP2over)[3:4]
  
## SEs of coefficients of binomial model recovered as:  
  {
    P2B <- rbind(c(1,0,-1,0),c(0,1,0,-1))
    (vcovP2B <- P2B \%*\% vcov(byP2over) \%*\% t(P2B))
  } 
  # practically equivalent to 
  vcov(byB)
  # and the original SEs:
  sqrt(diag(vcovP2B))

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
    list(blue ~ offset(.dynoffset) + 0, family = poisson())),
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

\dontrun{
#### Trinomial fit
(byP3 <- pois4mlogit(submodels = list(
    list(yellow ~ offset(.dynoffset) + 0, family = poisson()),
    list(blue ~ offset(.dynoffset) + 1, family = poisson()),
    list(purple ~ offset(.dynoffset) + 1, family = poisson())),
    data = toydata, types=c("yellow","blue", "purple")))

## Same fit with data in long form:
str(longdata <- reshape2long(toydata, types=c("yellow","purple","blue")))
(byP3l <- pois4mlogit(submodels = list(
    list(yellow ~ offset(.dynoffset) + 0, family = poisson()),
    list(blue ~ offset(.dynoffset) + 1, family = poisson()),
    list(purple ~ offset(.dynoffset) + 1, family = poisson())),
    data = longdata, types=c("yellow","blue", "purple")))

# Coefficients are identical but logL differ as explained for binomial example.
}
    
# Have a look at the 'See also' section for useful features 
# not illustrated in the present examples. 
% ____F I X M E____ illustrate 'aliases'?


}

