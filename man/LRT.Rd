\name{LRT}
\alias{LRT}
\title{
  Likelihood ratio tests of fixed and random effects.
}
\description{
\code{LRT} performs a likelihood ratio (LR) test between two model fits. It differs from another function with the same effect, \code{\link{fixedLRT}}, by its arguments (model fits for \code{LRT}, but all arguments required to fit the models for \code{fixedLRT}), and by the format of its return value. \code{LRT} determines which model is the more complete one by comparing model components including the fixed-effect, random-effect, residual-dispersion model specifications, and response families (offsets are ignored). Then, a standard test based on the asymptotic chi-square distribution is performed. In addition, parametric bootstrap p-values can be computed, either using the \emph{raw bootstrap} distribution of the likelihood ratio, or a bootstrap estimate of the Bartlett correction of the LR statistic. 

These different tests perform diffferently depending on the differences between the two models:\cr \cr 
 * If the models differ only by their fixed effects, the asymptotic LRT may be anticonservative, but the Bartlett-corrected one is generally well-calibrated. \cr \cr
 * If the two models differ by their random effects, tests based on the chi-square distribution (including their Bartlett-corrected version) may be poorly behaved, as such tests assume unbounded parameters, as contrasted to, e.g., necessarily positive variances. \cr
 In such cases the raw boostrap test may be the only reliable test. The procedure aims to detect and report such issues, but may not report all problems: users remain responsible for applying the tests in valid conditions (see \bold{Caveats} in Details section). In simple cases (such as comparing a fixed-effect to a mixed-effect model with the same fixed-effect term), the chi-square tests may not be reported. In other cases (see Examples) they may otherwise be reported, with a warning when the procedure detects some cases of estimates at the boundary for the full model, or detects cases where the LR statistic of bootstrap replicates is often zero (also suggesting that estimates are at the boundary in such replicates).\cr \cr
  * If the fits differ by the fixed effects terms of their residual-dispersion models (but not by any random effect specification), tests based on the chi-square distribution are reported. A bootstrap can be performed as in other cases.\cr
  * Tests for some cases of nested response families (e.g., the Poisson versus its extensions) are tentatively allowed.\cr
  * In some cases the full and the null models cannot be identified and the basic LRT based on the chi-square distribution will not be returned, but a bootstrap test may still be performed.\cr
  * The case where residual-dispersion models of either fit include random effects is problematic as, for such fits, the fitting procedure does not maximize the reported likelihood. The basic LRT is not returned when the two fits differ by their random effects, but is still performed otherwise (see Examples); and a bootstrap test may still be performed in both cases. 
}
\usage{
LRT(object, object2, boot.repl = 0L, resp_testfn = NULL, 
    simuland = eval_replicate, include="call", 
    #     many further arguments can be passed to spaMM_boot via the '...'
    #     These include arguments for parallel computations, such as
    # nb_cores, fit_env, 
    #     as well as other named arguments and spaMM_boot's own '...'
    ...)    
}
\arguments{
\item{object}{Fit object returned by a \pkg{spaMM} fitting function.}
\item{object2}{Optional second model fit to be compared to the first (their order does not matter, except in non-standard cases where the second model is taken as teh null one and a message is issued).}
\item{boot.repl}{the number of bootstrap replicates. }
\item{resp_testfn}{See argument \code{resp_testfn} of \code{\link{spaMM_boot}}.}
\item{simuland}{
  a function, passed to \code{\link{spaMM_boot}}. See argument \code{\link{eval_replicate}} for default value and requirements.
}
\item{include}{character string, either \code{"fit"} or  \code{"call"} (default; other values are ignored). 
Controls whether the fits compared, or only their calls, are included in the return value.}
%
\item{\dots}{
  Further arguments, passed to \code{\link{spaMM_boot}} (e.g., for parallelization).
}
}
\details{
  * \bold{Identifying a nested model}: \pkg{spaMM} tries to guess which is the nested model, but this is not always possible. 
  
  In particular, the fit with (substantially) lower logLik is not necessarily the one with fewer degrees of freedom, if the models are not nested, or even for some form of \dQuote{nestedness} with constrained parameters. The latter case may occur if one first fits a simple model, then fits a more complex model containing the first model \bold{but} one constrains estimates of shared parameters to the values given by the first fit: the second fit should then have a higher likelihood but it may still have fewer estimated parameters than the first one. A likelihood ratio test is invalid anyway when such constraints are used.  
  
  Second, although nestedness is relatively easy to assess for fixed effects, equivalent random-effect models can be specified by diverse syntaxes, so a simple textual comparison of the random-effect terms may not be enough, and model specifications that hinder such a comparison should be avoided. When differences in random effects are tested, the null distribution of the LR may include a probability mass in 1: the discussion in Details of \code{\link{get_RLRsim_args}} applies. 
  
  * \bold{\code{data} checks}: LRTs are invalid if the two models are fitted on different \bold{informative} data. Informative data may differ even if the input data are identical, because the fitting procedures ignore lines of data with missing information. The retained data may then differ for the two models, if additional variables of the more complete model are missing from some lines of data retained in the fit of the less complete model.  
  
  \code{LRT} perform some checks based on comparison of the data used by the two models. However, a perfect check all cases where the informative data differ is not easily defined. For example, data may appear different but still contain the same information if lines of the data are reordered. Hence, users remain ultimately responsible for checking whether informative data are equivalent in the two models. They should keep in mind that data processing is applied on a per-submodel basis for multivariate-response fits, and \code{\link{pois4mlogit}} has further criteria for definition of informative data for each model.

  * \bold{Bootstrap LRTs}: A raw bootstrap p-value can be computed from the simulated distribution as \code{(1+sum(t >= t0))/(N+1)} where \code{t0} is the original likelihood ratio, \code{t} the vector of bootstrap replicates and \code{N} its length. See Davison & Hinkley (1997, p. 141) for discussion of the adjustments in this formula. However, a computationally more economical use of the bootstrap is to provide a Bartlett correction for the likelihood ratio test in small samples. According to this correction, the mean value \eqn{m} of the likelihood ratio statistic under the null hypothesis is computed (here estimated by a parametric bootstrap) and the original LR statistic is multiplied by \eqn{n/m} where \eqn{n} is the number of degrees of freedom of the test. 
  
  * \bold{Caveats}: (1) An evaluated log-likelihood ratio can be slightly negative, e.g. when a fixed-effect model is compared to a mixed one, or a spatial random effect to a block effect, if parameters of the more complete model are estimated within bounds (e.g., variance>1e-06, or Matern smoothness>0.005) designed to avoid numerical singularities, while the less complete model corresponds to a boundary case (e.g., variance=0, or smoothness=0). The bootstrap procedure tries to identify these cases and then corrects slightly negative logL ratios to 0. (2) The Bartlett correction is applicable when the true distribution of the LRT departs smoothly from the chi-square distribution, but not in cases where it has a probability mass in zero (at typically occurs in the same boundary cases).     
}
\value{
\code{LRT} returns an object of class \code{fixedLRT}, unless some exception occurs.
This is a list with typical elements (depending on the options)
\item{fullfit, nullfit}{the fit objects for the full and null models, or only their calls, depending on the \code{include} argument;} 
\item{basicLRT}{A data frame including values of the likelihood ratio chi2 statistic, its degrees of freedom, and the p-value;}

and, if a bootstrap was performed: 
\item{rawBootLRT}{A data frame including values of the likelihood ratio chi2 statistic, its degrees of freedom, and the raw bootstrap p-value;}
\item{BartBootLRT}{A data frame including values of the Bartlett-corrected likelihood ratio chi2 statistic, its degrees of freedom, and its p-value;}
\item{bootInfo}{a list with the following elements:
 \describe{
\item{bootreps}{A table of fitted likelihoods for bootstrap replicates;  }
\item{meanbootLRT}{The mean likelihood ratio chi-square statistic for bootstrap replicates;  }
 } %- describe
} %- item

} %- value
\references{
 Bartlett, M. S. (1937) Properties of sufficiency and statistical tests. Proceedings of the Royal Society (London) A 160: 268-282.

 Davison A.C., Hinkley D.V. (1997) Bootstrap methods and their applications. Cambridge Univ. Press, Cambridge, UK.
}
\seealso{
See also \code{\link{fixedLRT}} and \code{\link{anova.HLfit}}\cr 
}
\examples{
## Using resp_testfn argument for bootstrap LRT:
\dontrun{
set.seed(1L)
d <- data.frame(success = rbinom(10, size = 1, prob = 0.9), x = 1:10)
xx <- cbind(1,d$x)
table(d$success)
m_x <- fitme(success ~ x, data = d, family = binomial())
m_0 <- fitme(success ~ 1, data = d, family = binomial())
#
# Bootstrap LRTs:
anova(m_x, m_0, boot.repl = 100,
      resp_testfn=function(y) {! is_separated(xx,as.numeric(y),verbose=FALSE)})
}

#### Various cases were asymptotic tests may be unreliable:

set.seed(123)
dat <- data.frame(g = rep(1:10, e = 10), x = (x<-rnorm(100)), 
                   y = 0.1 * x + rnorm(100))
m0 <- fitme(y ~ 1, data=dat) 

## (1) Models differing both by fixed and random effects: 
% so it's not the simple case where only the raw test is computed
#
# (note the warning for variance at boundary):
#
if (spaMM.getOption("example_maxtime")>11) { 
  m <- fitme(y ~ x + (1|g), data=dat)
  LRT(m,m0, boot.repl = 199L)
}
## See help("get_RLRsim_args") for a fast and accurate test procedure

## (2) Models differing also by residual-dispersion models:
#
if (spaMM.getOption("example_maxtime")>25) { 
  m <- fitme(y ~ x + (1|g), data=dat, resid.model= ~x)
  LRT(m,m0, boot.repl = 99L)
}

## (3) Models differing (also) by their random-effects in resid.model:
#
m <- fitme(y ~ x, data=dat, resid.model= ~1+(1|g)) 
LRT(m,m0) # no test performed


}
\keyword{ htest }
