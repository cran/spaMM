\name{LRT}
\alias{anova}
\alias{anova.HLfit}
\alias{LRT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Likelihood ratio test of fixed effects.
}
\description{
\code{LRT} performs a likelihood ratio (LR) test between two model fits, the ``full'' and the ``null'' model fits. The \code{anova} method for fit objects from spaMM is an alias for \code{LRT}. 

If the models differ only in their fixed effects, a standard test based on the asymptotic chi-square distribution is performed, with number of degrees of freedom determined by the function. In addition, parametric bootstrap p-values can be computed, either using the raw bootstrap distribution of the likelihood ratio, or a bootstrap estimate of the Bartlett correction of the LR statistic.
This function differs from \code{fixedLRT} in its arguments (model fits for \code{LRT}, but all arguments required to fit the models for \code{fixedLRT}), and in the format of its return value. By conceptual drift, \code{anova} works as an alias for \code{LRT}.  

If the two models differ in their random effects, a bootstrap test may be performed. Either the full and null model can be told apart by a simple comparison of the model formulas, or the model with the lower likelihood is assumed to be the null one (the latter is subject to numerical uncertainties when both fit are equivalent). No number of degrees of freedom is deduced, so no asymptotic test is performed. 

If the two models differ in both their fixed and random components, the same bootstrap test can be performed (see Examples), but the procedure further checks that the same model is nested in the other for both components. This requires that a simple comparison of the model formulas is sufficient to assess this for random effects, and may therefore fail). 
}
\usage{
\method{anova}{HLfit}(object, object2, ..., method="")
LRT(object, object2, boot.repl = 0, resp_testfn = NULL, 
    simuland = eval_replicate, 
    #     many further arguments can be passed to spaMM_boot via the '...'
    #     These include arguments for parallel computations, such as
    # nb_cores, fit_env,
    #     as well as other named arguments and spaMM_boot's own '...'
    ...)    
}
\arguments{
  \item{object,object2}{Two models fits being compared (their order does not matter).}
  \item{boot.repl}{the number of bootstrap replicates. }
\item{resp_testfn}{See argument \code{resp_testfn} of \code{\link{spaMM_boot}}}
  \item{method}{For development purposes, not documented. }

\item{simuland}{
  a function, passed to \code{\link{spaMM_boot}}. See argument \code{\link{eval_replicate}} for default value and requirements.
}  
\item{\dots}{
  Further arguments, passed to \code{\link{spaMM_boot}} (e.g., for parallelization). 
}
}
\details{
  A raw bootstrap p-value can be computed from the simulated distribution as \code{(1+sum(t >= t0))/(N+1)} where \code{t0} is the original likelihood ratio, \code{t} the vector of bootstrap replicates and \code{N} its length. See Davison & Hinkley (1997, p. 141) for discussion of the adjustments in this formula. However, a computationally more economical use of the bootstrap is to provide a Bartlett correction for the likelihood ratio test in small samples. According to this correction, the mean value \eqn{m} of the likelihood ratio statistic under the null hypothesis is computed (here estimated by a parametric bootstrap) and the original LR statistic is multiplied by \eqn{n/m} where \eqn{n} is the number of degrees of freedom of the test. 
  
  If random effects are tested, only the raw p-value is computed. Its null distribution may include a probability mass in 1 (the discussion in Details of \code{\link{get_RLRsim_args}} applies).
}
\value{
An object of class \code{fixedLRT}, actually a list with as-yet unstable format, but here with typical elements (depending on the options)
\item{fullfit}{the HLfit object for the full model;} 
\item{nullfit}{the HLfit object for the null model;}
\item{basicLRT}{A data frame including values of the likelihood ratio chi2 statistic, its degrees of freedom, and the p-value;}

and, if a bootstrap was performed: 
\item{rawBootLRT}{A data frame including values of the likelihood ratio chi2 statistic, its degrees of freedom, and the raw bootstrap p-value;}
\item{BartBootLRT}{A data frame including values of the Bartlett-corrected likelihood ratio chi2 statistic, its degrees of freedom, and its p-value;}
\item{bootInfo}{a list with the following elements:
 \describe{
\item{bootreps}{A table of fitted likelihoods for bootstrap replicates;  }
\item{meanbootLRT}{The mean likelihood ratio chi-square statistic for boostrap replicates;  }
 } %- describe
} %- item
} %- value
\references{
 Bartlett, M. S. (1937) Properties of sufficiency and statistical tests. Proceedings of the Royal Society (London) A 160: 268-282.

 Davison A.C., Hinkley D.V. (1997) Bootstrap methods and their applications. Cambridge Univ. Press, Cambridge, UK.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
See also \code{\link{fixedLRT}} and \code{\link{get_RLRsim_args}}.
}
\examples{
data("wafers")
## Gamma GLMM with log link
m1 <- HLfit(y ~X1+X2+X1*X3+X2*X3+I(X2^2)+(1|batch),family=Gamma(log),
          resid.model = ~ X3+I(X3^2) ,data=wafers,method="ML")
m2 <- update(m1,formula.= ~ . -I(X2^2))
anova(m1,m2)

# Using resp_testfn argument:
\dontrun{
set.seed(1L)
d <- data.frame(success = rbinom(10, size = 1, prob = 0.9), x = 1:10)
xx <- cbind(1,d$x)
table(d$success)
m_x <- fitme(success ~ x, data = d, family = binomial())
m_0 <- fitme(success ~ 1, data = d, family = binomial())
anova(m_x, m_0, boot.repl = 100,
      resp_testfn=function(y) {! is_separated(xx,as.numeric(y),verbose=FALSE)})
}

## Models differing both in fixed and random effects:
if (spaMM.getOption("example_maxtime")>11) { 
 set.seed(123)
 dat <- data.frame(g = rep(1:10, e = 10), x = (x<-rnorm(100)), 
                   y = 0.1 * x + rnorm(100))
 m <- fitme(y ~ x + (1|g), data=dat)
 m0 <- fitme(y ~ 1, data=dat) 
 (bootpval <- LRT(m,m0, boot.repl = 199L)$rawBootLRT$p_value)
 ## See help("get_RLRsim_args") for a fast and accurate test procedure
}

}
\keyword{ htest }
