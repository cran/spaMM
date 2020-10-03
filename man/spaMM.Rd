\name{spaMM}
% ----- \Rdversion{1.1}
\alias{spaMM}
\alias{spaMM-package}

\docType{package}

\title{Inference in mixed models, in particular spatial GLMMs}

\description{Fits a range of mixed-effect models, including those with spatially correlated random effects. The random effects are either Gaussian (which defines GLMMs), or other distributions (which defines the wider class of hierarchical GLMs), or simply absent (which makes a LM or GLM).  
spaMM is designed to be used through the high-level fitting functions 
\code{\link{fitme}} (the most general function), \code{\link{HLfit}} (sometimes faster, for non-spatial models), \code{\link{HLCor}} (sometimes faster, for conditional-autoregressive models and fixed-correlation models), \code{\link{corrHLfit}} (now of lesser interest); and additional functions such as \code{\link{fixedLRT}} for likelihood-ratio testing, \code{\link{simulate}} and \code{\link{predict}}.

Both maximum likelihood (ML) and restricted likelihood (REML) can be used for linear mixed models, and extensions of these methods using Laplace approximations are used for non-Gaussian random response. Several variants of these methods discussed in the literature are included (see Details in \code{\link{HLfit}}), the most notable of which may be \dQuote{PQL/L} for binary-response GLMMs (see Example for \code{\link{arabidopsis}} data). PQL methods implemented in spaMM are closer to (RE)ML methods than those implemented in \code{MASS::glmmPQL}. 
}

\details{
  The standard response families \code{gaussian}, \code{binomial}, \code{poisson}, and \code{Gamma} are handled, as well as negative binomial (see \code{\link{negbin}}), zero-truncated poisson and negative binomial, and Conway-Maxwell-Poisson response (see \code{\link{Tpoisson}}, \code{\link{Tnegbin}} and \code{\link{COMPoisson}}). A \code{multi} family look-alike is also available for \code{\link{multinomial}} response, with some constraints. 
  
  The variance parameter of residual error is denoted \eqn{\phi} (\code{phi}): this is the residual variance for gaussian response, but for Gamma-distributed response, the residual variance is \eqn{\phi}\eqn{\mu^2} where \eqn{\mu} is expected response. A fixed-effects linear predictor for \eqn{\phi}, modeling heteroscedasticity, can be considered (see Examples). 
  
  The package fits models including several nested or crossed random effects, including autocorrelated ones with the following correlation models: \code{\link{Matern}}, \code{\link{Cauchy}}, interpolated Markov Random Fields (\code{\link{IMRF}}), first-order autoregressive (\code{\link{AR1}}), conditional autoregressive as specified by an \code{\link{adjacency}} matrix, or any fixed correlation matrix (\code{\link{corrMatrix}}). GLMMs and HGLMs are fit via Laplace approximations for (1) the marginal likelihood with respect to random effects and (2) the restricted likelihood (as in REML), i.e. the likelihood of random effect parameters given the fixed effect estimates. 
  
     All handled models can be formulated in terms of a linear predictor of the traditional form \code{offset}+ \bold{X}\eqn{\beta} + \bold{Z b}, where \bold{X} is the design matrix of fixed effects, \eqn{\beta} (\code{beta}) is a vector of fixed-effect coefficients, \bold{Z} is a \dQuote{design matrix} for the random effects (which is instead denoted \bold{M}=\bold{ZAL} elsewhere in the package documentation), and \bold{b} a vector of random effect values. The general structure of \code{Mb} is described in \code{\link{random-effects}}.
  
   The syntax for formulas extends that used in the \code{lme4} package. In particular, \bold{non-autocorrelated} random effects are specified using the \code{(1|<block>)} syntax, and \emph{Gaussian} random-coefficient terms by \code{(<rhs>|<block>)}.The double-vertical syntax, \code{( rhs || lhs )}, is interpreted as in \code{lme4}. Any such term is immediately converted to \code{( ( 1 | lhs ) + ( 0 + lhs | rhs ) )}, and should be counted as two random effects for all purposes (e.g., for fixing the variances of the random effects). \bold{Autocorrelated} random effects are specified by adding some prefix to this syntax, \code{<prefix>(1|.)}, e.g. \code{Matern(1|long+lat)}. 
   
   Since version 2.6.0, it is possible to fit some \dQuote{some autocorrelated random-coefficient} models by a syntax consistent with that of random-coefficient terms, \code{<prefix>(<rhs>|.)}. For example, independent Matérn effects can be fitted for males and females by using the syntax \code{Matern(male|.) + Matern(female|.)}, where \code{male} and \code{female} are TRUE/FALSE factors. A numeric variable \code{z} can also be considered, in which case the proper syntax is \code{Matern(0+z|.)}, which represents an autocorrelated random-slope (only) term (or, equivalently, a direct specification of heteroscedasticy of the Matérn random effect). All these effects are achieved by direct control of the elements of the incidence matrix \bold{Z} through the \code{<rhs>} term: for numeric \code{z}, such elements are multiplied by \code{z} values, and thus provide a variance of order O(\code{z} \bold{squared}). By contrast, \code{Matern(z|.)} is not defined. It could mean that a correlation structure between random intercepts and random slopes is to be combined with a Matérn correlation structure, but no way of combining them is yet defined and implemented.  
% cf fns .calc_Z_model_matrix() and .calc_raw_ZA().
% tests: test-spaMM.R and in tests_private/test-syntax-Matern-RHS.R

Since version 2.7.0, the syntax \code{(z-1|.)}, for \bold{numeric} \code{z} only, can also be used to fit \bold{some heteroscedastic non-Gaussian} random effects. For example, a Gamma random-effect term \code{(wei-1|block)} specifies an heteroscedastic Gamma random effect \eqn{u} with constant mean 1 and variance \code{wei^2} \eqn{\lambda}, where \eqn{\lambda} is still the estimated variance parameter. See Details of \code{\link{negbin}} for a possible application. Here, this effect is not implemented through direct control of \bold{Z} (multiplying the elements of an incidence matrix \bold{Z} by \code{wei}), as this would have a different effect on the distribution of the random effect term. \code{(z|.)} is not defined for \emph{non-Gaussian} random effects. It could mean that a correlation structure between random intercepts and random slopes for (say) Gamma-distributed random effects is considered, but such correlation structures are not well-specified by their correlation matrix. 

}
\author{
\code{spaMM} was initially published by \enc{François}{Francois} Rousset and Jean-Baptiste Ferdy, and is continually developed by F. Rousset and tested by Alexandre Courtiol.  
}
\references{
Lee, Y., Nelder, J. A. and Pawitan, Y. (2006). Generalized linear models with random effects: unified analysis via
h-likelihood. Chapman & Hall: London.

Rousset F., Ferdy, J.-B. (2014) Testing environmental and genetic effects in the presence of spatial autocorrelation. Ecography, 37: 781-790.
    \url{http://dx.doi.org/10.1111/ecog.00566}

}
\seealso{
The \code{test} directory of the package provides many additional examples of spaMM usage beyond those from the formal documentation. 
}

\examples{
data("wafers")
data("scotlip") ## loads 'scotlip' data frame, but also 'Nmatrix'

##     Linear model
fitme(y ~ X1, data=wafers)

##     GLM
fitme(y ~ X1, family=Gamma(log), data=wafers)
fitme(cases ~ I(log(population)), data=scotlip, family=poisson)

##     Non-spatial GLMMs
fitme(y ~ 1+(1|batch), family=Gamma(log), data=wafers)
fitme(cases ~ 1+(1|gridcode), data=scotlip, family=poisson)
#
# Random-slope model (mind the output!)        
fitme(y~X1+(X2|batch),data=wafers, method="REML")

## Spatial, conditional-autoregressive GLMM
if (spaMM.getOption("example_maxtime")>2) {   
  fitme(cases ~ I(log(population))+adjacency(1|gridcode), data=scotlip, family=poisson, 
        adjMatrix=Nmatrix) # with adjacency matrix provided by data("scotlip")
} 
# see ?adjacency for more details on these models 

## Spatial, geostatistical GLMM: 
# see e.g. examples in ?fitme, ?corrHLfit, ?Loaloa, or ?arabidopsis;
# see examples in ?Matern for group-specific spatial effects.

##     Hierachical GLMs with non-gaussian random effects
 data("salamander")
if (spaMM.getOption("example_maxtime")>1) {   
 # both gaussian and non-gaussian random effects
 fitme(cbind(Mate,1-Mate)~1+(1|Female)+(1|Male),family=binomial(),
        rand.family=list(gaussian(),Beta(logit)),data=salamander)
 
 # Random effect of Male nested in that of Female:
 fitme(cbind(Mate,1-Mate)~1+(1|Female/Male),
       family=binomial(),rand.family=Beta(logit),data=salamander)
 # [ also allowed is cbind(Mate,1-Mate)~1+(1|Female)+(1|Male \%in\% Female) ]
}

##    Modelling residual variance ( = structured-dispersion models)    
# GLM response, fixed effects for residual variance 
fitme( y ~ 1,family=Gamma(log),
      resid.model = ~ X3+I(X3^2) ,data=wafers)
#
# GLMM response, and mixed effects for residual variance
if (spaMM.getOption("example_maxtime")>1.5) {   
  fitme(y ~ 1+(1|batch),family=Gamma(log),
        resid.model = ~ 1+(1|batch) ,data=wafers)
}
}
\keyword{ package }
\keyword{spatial}
% should add \keyword{smooth} later... and htest ? models ? see RShowDoc("KEYWORDS")
