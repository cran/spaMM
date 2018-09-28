\name{AIC}
\alias{get_any_IC}
\alias{AIC}
\alias{AIC.HLfit}
\alias{extractAIC}
\alias{extractAIC.HLfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extractors for informations criteria such as AIC
}
\description{
\code{get_any_IC} computes model selection/information criteria such as AIC. See Details for more information about these criteria. The other extractors \code{AIC} and \code{extractAIC} are methods for \code{HLfit} objects of generic functions defined in other packages; \code{AIC} is equivalent to \code{get_any_IC}, and \code{extractAIC} additionally returns a number of degrees of freedom. 
}
\usage{
get_any_IC(object, ..., verbose=interactive() ,also_cAIC=TRUE)
\method{AIC}{HLfit}(object, ..., k, verbose=interactive() ,also_cAIC=TRUE)
\method{extractAIC}{HLfit}(fit, scale, k, ..., verbose=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object, fit}{A object of class \code{HLfit}, as returned by the fitting functions in \code{spaMM}.}
  \item{scale, k}{Currently ignored, but are required in the definitions for consistency with the generic.}
  \item{verbose}{ Whether to print the model selection criteria or not. }
  \item{also_cAIC}{Whether to include the conditional AIC in the result (its computation may be slow).}
  \item{\dots}{Other arguments that may be needed by some method.}
}
\details{\code{get_any_IC} computes, optionally prints, and returns invisibly the following quantities. The \bold{conditional AIC} (Vaida and Blanchard 2005) is a relative measure of quality of prediction of new realizations of a mixed model, conditional on the realized values of the random effects. It involves the conditional likelihood, and degrees of freedom for (i) estimated residual error parameters and (ii) the overall linear predictor characterized by the \bold{Effective degrees of freedom} already discussed by previous authors including Lee and Nelder (2001), which gave a general formula for it in HGLMs. 
Both a plug-in \dQuote{asymptotic} estimate of the conditional AIC and of this effective df are returned by \code{get_any_IC}.  Note that these may be biased estimates of conditional AIC and effective df, and that more refined formulas are discussed in the literature (e.g. Overholser and Xu 2014), some of which may be implemented in future versions of \code{get_any_IC}. 
Lee et al. (2006) and Ha et al (2007) defined a corrected AIC [i.e., AIC(D*) in their eq. 7] which is here interpreted as the conditional AIC. Also returned are the \bold{marginal AIC} (Akaike's classical AIC), and a focussed AIC for dispersion parameters (\bold{dispersion AIC}) discussed by Ha et al (2007; eq.10). This diversity of criteria should encourage users to think twice before applying model selection automatically, which is no better although more fashionable than misuses of simple null hypothesis testing. Also, alternative procedures for model choice can be considered (e.g. Cox and Donnelly, 2011, p. 130-131).
}
\value{
  For \code{AIC} and \code{get_any_IC}, a numeric vector whose elements are described in the Details.
  
  For \code{extractAIC}, a numeric vector of length 2, with first and second elements giving

  \item{edf}{the degree of freedom of the fixed-effect terms of the model
    for the fitted model \code{fit}.}

\item{AIC}{the (marginal) Akaike Information Criterion for \code{fit}.}
}
\references{
Cox, D. R. and Donnelly C. A. (2011) Principles of Applied Statistics. Cambridge Univ. Press.

Ha, I. D., Lee, Y. and MacKenzie, G. (2007) Model selection for multi-component frailty models. Statistics in Medicine 26: 4790-4807.

Overholser R., and Xu R. (2104) Effective degrees of freedom and its application to conditional AIC for linear mixed-effects models with correlated error structures. J. Multivariate Anal. 132: 160-170.

Vaida, F., and Blanchard, S. (2005) Conditional Akaike information for mixed-effects models. Biometrika 92, 351-370.
}
\examples{
data("wafers")
m1 <- HLfit(y ~X1+X2+(1|batch),
          resid.model = ~ 1 ,data=wafers,HLmethod="ML")
get_any_IC(m1) 
extractAIC(m1)
}
\keyword{models}
