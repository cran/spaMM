\name{LRT}
\alias{anova}
\alias{anova.HLfit}
\alias{anova.LMLT}
\alias{LRT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  ANOVA tables, and likelihood ratio tests of fixed and random effects.
}
\description{
The \code{anova} method for fit objects from \pkg{spaMM} has two uses: if a single fit object is provided, ANOVA tables may be returned, with specific procedures for univariate-response LMs, GLMs and LMMs (see Details). Alternatively, if a second fit object is provided (\code{object2} argument), \code{anova} performs as an alias for \code{LRT}.  The \code{LRT} function here differs from \code{\link{fixedLRT}} by its arguments (model fits for \code{LRT}, but all arguments required to fit the models for \code{fixedLRT}), and by the format of its return value.

\code{LRT} performs a likelihood ratio (LR) test between two model fits, the ``full'' and the ``null'' model fits. It determines which model is the more complete one by comparing model components including the fixed-effect, random-effect, residual-dispersion model specifications, and response families (offsets are ignored). Then, a standard test based on the asymptotic chi-square distribution is performed. In addition, parametric bootstrap p-values can be computed, either using the \emph{raw bootstrap} distribution of the likelihood ratio, or a bootstrap estimate of the Bartlett correction of the LR statistic. 

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
\method{anova}{HLfit}(object, object2, type = "2", method="", ...)
#
LRT(object, object2, boot.repl = 0L, resp_testfn = NULL, 
    simuland = eval_replicate, 
    #     many further arguments can be passed to spaMM_boot via the '...'
    #     These include arguments for parallel computations, such as
    # nb_cores, fit_env,
    #     as well as other named arguments and spaMM_boot's own '...'
    ...)    
}
\arguments{
\item{object}{Fit object returned by a \pkg{spaMM} fitting function.}
\item{object2}{Optional second model fit to be compared to the first (their order does not matter, except in non-standard cases where the second model is taken as teh null one and a message is issued).}
\item{type}{ANOVA type for LMMs, as interpreted by \pkg{lmerTest}. Note that the default (single-term deletion ANOVA) differs from that of \pkg{lmerTest}.} 
%
\item{boot.repl}{the number of bootstrap replicates. }
\item{resp_testfn}{See argument \code{resp_testfn} of \code{\link{spaMM_boot}}.}
\item{simuland}{
  a function, passed to \code{\link{spaMM_boot}}. See argument \code{\link{eval_replicate}} for default value and requirements.
}
\item{method}{
  Only non-default value is \code{"t.Chisq"} which forces evaluation of a table of chi-squared tests for each fixed-effect term, using the classical \dQuote{Wald} test (see Details). Further methods are available through the \dots for specific classes of models.
} 
%
\item{\dots}{
  Further arguments, passed to \code{\link{spaMM_boot}} (e.g., for parallelization) in the case of LRTs. For ANOVA tables, arguments of functions \code{\link[stats]{anova.lm}} \code{\link[stats]{anova.glm}}, 
  and \code{\link{as_LMLT}}, % including, notably, its \code{check} argument
  respectively for LMs, GLMs and LMMs, may be handled (e.g. the \code{test} argument for \code{anova.glm}).
}
}
\details{
  *   The \bold{ANOVA-table} functionality has been included here mainly to provide access to F tests (including, for LMMs, the \dQuote{Satterthwaite method} as developed by Fai and Cornelius, 1996), using pre-existing procedures as template or backend for expediency and familiarity:\cr 
  \enumerate{
 \item ANOVA tables \bold{for LMs and GLMs} have been conceived to replicate the functionality, output format and details of base R \code{anova}, and therefore replicate some of their limitations, e.g., they only perform sequential analysis (\dQuote{type 1}) in the same way as \code{anova.lm} and \code{anova.glm}. However, a difference occurs for Gamma GLMs, because the dispersion estimates for Gamma GLMs differ between \code{stats::glm} and \pkg{spaMM} fits (see Details in \code{\link{method}}). Therefore, F tests and Mallows' Cp differ too; results from \pkg{spaMM} REML fits being closer than ML fits to those from \code{glm()} fits;\cr \item \bold{For LMMs}, ANOVA tables are provided by interfacing \code{lmerTest::anova} (with non-default \code{type}). This procedure should handle all types of LMMs that can be fitted by \pkg{spaMM}; yet, the displayed information should be inspected to check that some fitted random-effect parameters are not ignored when computing information for the Satterthwaite method. 
\item  For fitted models that do not lay within previous categories, such as \bold{GLMMs}, models with a \bold{residual-dispersion} submodel, and \bold{multivariate-response} models, a table of tests for single-term deletions using the classical \dQuote{Wald} chi-squared test based on coefficient values and their conditional standard error estimates will be returned. LRTs (moreover, with bootstrap correction) are more reliable than such tests and, as calling them requires a second model to be explicitly specified, they may also help users thinking about the hypothesis they are testing. 
  }

  * \bold{Bootstrap LRTs}: A raw bootstrap p-value can be computed from the simulated distribution as \code{(1+sum(t >= t0))/(N+1)} where \code{t0} is the original likelihood ratio, \code{t} the vector of bootstrap replicates and \code{N} its length. See Davison & Hinkley (1997, p. 141) for discussion of the adjustments in this formula. However, a computationally more economical use of the bootstrap is to provide a Bartlett correction for the likelihood ratio test in small samples. According to this correction, the mean value \eqn{m} of the likelihood ratio statistic under the null hypothesis is computed (here estimated by a parametric bootstrap) and the original LR statistic is multiplied by \eqn{n/m} where \eqn{n} is the number of degrees of freedom of the test. 
  
  When models differ by their random-effect specifications, distinguishing the full from the null model is not easy. In particular, equivalent models can be specified by diverse syntaxes, so a simple textual comparison of the random-effect terms may not be enough, and model specifications that hinder such a comparison should be avoided. When differences in random effects are tested, the null distribution of the LR may include a probability mass in 1: the discussion in Details of \code{\link{get_RLRsim_args}} applies. 
  
  * \bold{Caveats}: (1) An evaluated log-likelihood ratio can be slightly negative, e.g. when a fixed-effect model is compared to a mixed one, or a spatial random effect to a block effect, if parameters of the more complete model are estimated within bounds (e.g., variance>1e-06, or Matern smoothness>0.005) designed to avoid numerical singularities, while the less complete model corresponds to a boundary case (e.g., variance=0, or smoothness=0). The bootstrap procedure tries to identify these cases and then corrects slightly negative logL ratios to 0. (2) The Bartlett correction is applicable when the true distribution of the LRT departs smoothly from the chi-square distribution, but not in cases where it has a probability mass in zero (at typically occurs in the same boundary cases).     
}
\value{
\code{LRT} returns an object of class \code{fixedLRT}, actually a list with typical elements (depending on the options)
\item{fullfit}{the HLfit object for the full model;} 
\item{nullfit}{the HLfit object for the null model;}
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

When ANOVA tables are computed, the return format is that of the function called (\code{lmerTest::anova} for LMMs) or emulated (for LMs or GLMs). For GLMs, by default no test is reported, as has been the default for \code{anova.glm} before R 4.4.0.
} %- value
\references{
 Bartlett, M. S. (1937) Properties of sufficiency and statistical tests. Proceedings of the Royal Society (London) A 160: 268-282.

 Davison A.C., Hinkley D.V. (1997) Bootstrap methods and their applications. Cambridge Univ. Press, Cambridge, UK.
 
Fai AH, Cornelius PL (1996). Approximate F-tests of multiple degree of freedom hypotheses in generalised least squares analyses of unbalanced split-plot experiments.
Journal of Statistical Computation and Simulation, 54(4), 363-378. \doi{10.1080/00949659608811740}
}
\seealso{
See also \code{\link{fixedLRT}} for a different interface to LRTs,\cr 
\code{\link{get_RLRsim_args}} for efficient simulation-based implementation of exact likelihood ratio tests for testing the presence of variance components,\cr
\code{\link{as_LMLT}} for the interface to \code{lmerTest::anova},\cr
and \code{\link{summary.HLfit}}\code{(.,details=list(<true|"Wald">))} for reporting the p-value for each t-statistic in the summary table for fixed effects, either by Student's t distribution, or by the approximation of t^2 distribution by the Chi-squared distribution (\dQuote{Wald's test}).
}
\examples{
data("wafers")
## Gamma GLMM with log link
m1 <- HLfit(y ~X1+X2+X1*X3+X2*X3+I(X2^2)+(1|batch),family=Gamma(log),
          resid.model = ~ X3+I(X3^2) ,data=wafers,method="ML")
m2 <- update(m1,formula.= ~ . -I(X2^2))
#
anova(m1,m2) # LRT

## 'anova' (Wald chi-squared tests...) for GLMM or model with a 'resid.model'
anova(m1) 

## ANOVA table for GLM
# Gamma example, from McCullagh & Nelder (1989, pp. 300-2), as in 'glm' doc:
clotting <- data.frame(
    u = c(5,10,15,20,30,40,60,80,100),
    lot1 = c(118,58,42,35,27,25,21,19,18),
    lot2 = c(69,35,26,21,18,16,13,12,12))
spglm <- fitme(lot1 ~ log(u), data = clotting, family = Gamma, method="REML")
anova(spglm, test = "F") 
anova(spglm, test = "Cp") 
anova(spglm, test = "Chisq")
anova(spglm, test = "Rao") 

## ANOVA table for LMM
if(requireNamespace("lmerTest", quietly=TRUE)) {
  lmmfit <- fitme(y ~X1+X2+X1*X3+X2*X3+I(X2^2)+(1|batch),data=wafers)
  print(anova(lmmfit)) # => Satterthwaite method, here giving p-values 
                       #   quite close to traditional t-tests given by:
  summary(lmmfit, details=list(p_value=TRUE))
}

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
