\name{AIC}
\alias{get_any_IC}
\alias{AIC}
\alias{AIC.HLfit}
\alias{extractAIC}
\alias{extractAIC.HLfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extractors for information criteria such as AIC
}
\description{
\code{get_any_IC} computes model selection/information criteria such as AIC. See Details for more information about these criteria. The other extractors \code{AIC} and \code{extractAIC} are methods for \code{HLfit} objects of generic functions defined in other packages: \code{AIC} is equivalent to \code{get_any_IC}, and \code{extractAIC} returns the marginal AIC and the number of degrees of freedom for the fixed effects. 
}
\usage{
get_any_IC(object, nsim=0L, ..., verbose=interactive(),
           also_cAIC=TRUE, short.names=NULL)
\method{AIC}{HLfit}(object, nsim=0L, ..., k, verbose=interactive(),
                    also_cAIC=TRUE, short.names=NULL)
\method{extractAIC}{HLfit}(fit, scale, k, ..., verbose=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object, fit}{A object of class \code{HLfit}, as returned by the fitting functions in \code{spaMM}.}
  \item{scale, k}{Currently ignored, but are required in the definitions for consistency with the generic.}
  \item{verbose}{ Whether to print the model selection criteria or not. }
  \item{also_cAIC}{Whether to include the plug-in estimate of conditional AIC in the result (its computation may be slow).}
  \item{nsim}{Controls whether to include the bootstrap estimate of conditional AIC (see Details) in the result. If positive, \code{nsim} gives the number of bootstrap replicates.}
  \item{short.names}{NULL, or boolean; controls whether the return value uses short names (\code{mAIC}, etc., as shown by screen output if \code{verbose} is TRUE), or the descriptive names (\code{"   marginal AIC:"}, etc.) also shown in the screen output. Short names are more appropriate for programming but descriptive names may be needed for back-compatibility. The default (NULL) ensures back-compatibility by using descriptive names unless the bootstrap estimate of conditional AIC is reported.}
  \item{\dots}{Other arguments that may be needed by some method. If \code{nsim} is positive, a \code{seed} argument may be passed to \code{simulate}, and the other \dQuote{\dots} will be used in a call to \code{\link{dopar}} to control the optional parallel execution of the bootstrap computations.
  }
}
\details{
The AIC is a measure (by Kullback-Leibler directed distance, up to an additive constant) of quality of prediction of new data by a fitted model. 
Comparing information criteria may be viewed as a fast alternative to a comparison of the predictive accuracy of different models by cross-validation. Further procedures for model choice may also be useful (e.g. Williams, 1970; Lewis et al. 2010).

The \bold{conditional AIC} (Vaida and Blanchard 2005) applies the AIC concept to new realizations of a mixed model, conditional on the realized values of the random effects. Lee et al. (2006) and Ha et al (2007) defined a corrected AIC [i.e., AIC(D*) in their eq. 7] which is here interpreted as the conditional AIC. 

Such Kullback-Leibler relative distances cannot generally be evaluated exactly and various estimates have been discussed.
\code{get_any_IC} computes, optionally prints, and returns invisibly one or more of the following quantities: (1) Akaike's classical AIC (\bold{marginal AIC}, \code{mAIC}); (2) a plug-in estimate (\code{cAIC}) and/or a bootstrap estimate (\code{b_cAIC}) of the conditional AIC; and (3) a focussed AIC for dispersion parameters (\bold{dispersion AIC}, \code{dAIC}). 

For the \bold{conditional AIC}, Vaida and Blanchard's plug-in estimator involves the conditional likelihood, and degrees of freedom for (i) estimated residual error parameters and (ii) the overall linear predictor characterized by the \bold{Effective degrees of freedom} already discussed by previous authors including Lee and Nelder (1996), which gave a plug-in estimator (\eqn{p_D}) for it in HGLMs. 
By default, the plug-in estimate of both the conditional AIC and of \eqn{n-p_D} (\code{GoFdf}, where \eqn{n} is the length of the response vector) are returned by \code{get_any_IC}. But these are biased estimates of conditional AIC and effective df, and an alternative procedure is available if a non-default positive \code{nsim} value is used. In that case, the conditional AIC is estimated by a bootstrap version of Saefken et al. (2014)'s equation 2.5; this involves refitting the model to each bootstrap samples, so it may take time, and a full cross-validation procedure might as well be considered for model selection. 

The dispersion AIC has been defined from restricted likelihood by Ha et al (2007; eq.10). The present implementation will use restricted likelihood only if made available by an REML fit, otherwise marginal likelihood is used.
}
\value{
  For \code{AIC} and \code{get_any_IC}, a numeric vector whose possible elements are described in the Details, and whose names are controlled by the \code{short.names} argument. Note that the bootstrap computation actually makes sense and works also for fixed-effect models (although it is not clear how useful it is in that case). The return value will still refer to its results as conditional AIC.
  
  For \code{extractAIC}, a numeric vector of length 2, with first and second elements giving

  \item{edf}{the degree of freedom of the fixed-effect terms of the model
    for the fitted model \code{fit}.}

  \item{AIC}{the (marginal) Akaike Information Criterion for \code{fit}.}

This output aims to be equivalent (except for the explicit names) to the one from \code{stats::extractAIC}, despite the obscurities of the latter's documentation, and is indeed equivalent in value for GLMs (see Examples). 
}
\references{
%Cox, D. R. and Donnelly C. A. (2011) Principles of Applied Statistics. Cambridge Univ. Press.

Ha, I. D., Lee, Y. and MacKenzie, G. (2007) Model selection for multi-component frailty models. Statistics in Medicine 26: 4790-4807.

Lee Y. and Nelder. J. A. 1996. Hierarchical generalized linear models (with discussion). J. R. Statist. Soc. B, 58: 619-678. 

Lewis, F., Butler, A. and Gilbert, L. (2011), A unified approach to model selection using the likelihood ratio test. Methods in Ecology and Evolution, 2: 155-162. \doi{10.1111/j.2041-210X.2010.00063.x}

%Overholser R., and Xu R. (2104) Effective degrees of freedom and its application to conditional AIC for linear mixed-effects models with correlated error structures. J. Multivariate Anal. 132: 160-170.

Saefken B., Kneib T., van Waveren C.-S., Greven S. (2014) A unifying approach to the estimation of the conditional Akaike information in generalized linear mixed models. Electron. J. Statist. 8, 201-225. 

Vaida, F., and Blanchard, S. (2005) Conditional Akaike information for mixed-effects models. Biometrika 92, 351-370.

Williams D.A. (1970) Discrimination between regression models to determine the pattern of enzyme synthesis in synchronous cell cultures. Biometrics 26: 23-32.
}
\examples{
data("wafers")
m1 <- fitme(y ~ X1+X2+X3+X1*X3+X2*X3+I(X2^2)+(1|batch), data=wafers, 
            family=Gamma(log))

get_any_IC(m1) 
# => The plug-in estimate is stored in the 'm1' object 
#    as a result of the previous computation, and is now returned even by: 
get_any_IC(m1, also_cAIC=FALSE)

if (spaMM.getOption("example_maxtime")>4) {
 get_any_IC(m1, nsim=100L, seed=123) # provides bootstrap estimate of cAIC.
 # (parallelisation options could be used, e.g. nb_cores=detectCores()-1L)
}

extractAIC(m1)

\dontrun{
# Checking consistency with glm example from help("stats::extractAIC"):
utils::example(glm)
dataf <- data.frame(counts=counts,outcome=outcome, treatment=treatment)
extractAIC(glm.D93) 
extractAIC(fitme(counts ~ outcome + treatment, family = poisson(), data=dataf))
lm.D93 <- glm(counts ~ outcome + treatment)
extractAIC(lm.D93)
extractAIC(fitme(counts ~ outcome + treatment, data=dataf))
}

}
\keyword{models}
