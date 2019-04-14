\name{spaMM}
% ----- \Rdversion{1.1}
\alias{spaMM}
\alias{spaMM-package}

\docType{package}

\title{Inference in mixed models, in particular spatial GLMMs}

\description{Fits a range of mixed-effect models, including those with spatially correlated random effects. The random effects are either Gaussian (which defines GLMMs), or other distributions (which defines the wider class of hierarchical GLMs), or simply absent (which makes a GLM).  
}

\details{

  The standard response families \code{gaussian}, \code{binomial}, \code{poisson}, and \code{Gamma} are handled, as well as negative binomial (see \code{\link{negbin}}), zero-truncated poisson and negative binomial, and Conway-Maxwell-Poisson response (see \code{\link{Tpoisson}}, \code{\link{Tnegbin}} and \code{\link{COMPoisson}}). A \code{multi} family look-alike is also available for \code{\link{multinomial}} response, with some constraints. The variance parameter of residual error is denoted \eqn{\phi} (\code{phi}): this is the residual variance for gaussian response, but for Gamma-distributed response, the residual variance is \eqn{\phi}\eqn{\mu^2} where \eqn{\mu} is expected response. A fixed-effects linear predictor for \eqn{\phi}, modeling heteroscedasticity, can be considered (see Examples). 
  
  The package fits models including several nested or crossed random effects, including autocorrelated ones with the following correlation models: \code{\link{Matern}}, \code{\link{Cauchy}}, interpolated Markov Random Fields (\code{\link{IMRF}}), first-order autoregressive (\code{\link{AR1}}), conditional autoregressive as specified by an \code{\link{adjacency}} matrix, or any fixed correlation matrix (\code{\link{corrMatrix}}). GLMMs and HGLMs are fit via Laplace approximations for (1) the marginal likelihood with respect to random effects and (2) the restricted likelihood (as in REML), i.e. the likelihood of random effect parameters given the fixed effect estimates. 
  
     All handled models can be formulated in terms of a linear predictor of the standard form \code{offset}+ \bold{X}\eqn{\beta} + \bold{Z b}, where \bold{X} is the design matrix of fixed effects, \eqn{\beta} (\code{beta}) is a vector of fixed-effect coefficients, \bold{Z} is a design matrix of random effects (typically an incidence matrix), and \bold{b} a vector of random effect values. 
  
  The structure of the random effects \code{Zb} can generally be described by the following steps. First, independent and  identically distributed (iid) random effects \bold{u} are drawn from one of the following distributions: gaussian, Beta-distributed, Gamma and inverse-Gamma distributed random effects, implemented as detailed in the \code{\link{HLfit}} documentation. The variance(s) of random effects (\eqn{u}) is (are) denoted \eqn{\lambda}. Second, a transformation \bold{v}\eqn{=f}(\bold{u}) is applied (\bold{v} elements are still iid). Third, correlated random effects are obtained as \bold{Mv}, where the matrix \bold{M} can describe spatial correlation between observed locations, block effects (or repeated observations in given locations), and correlations involving unobserved locations. See Details in \code{\link{covStruct}} for the general form of \bold{M} as a matrix product \bold{ZAL}. In most cases \bold{M} is determined from the model formula, but it can also be input directly (e.g., to describe genetic correlations).     
  
   The syntax for formulas extends that used in the \code{lme4} package. In particular, non-autocorrelated random effects are specified using the \code{(1|<block>)} syntax, and \emph{Gaussian} random-coefficient terms by \code{(<rhs>|<block>)}. Autocorrelated random effects are specified by adding some prefix to this syntax, \code{<prefix>(1|.)}, e.g. \code{Matern(1|long+lat)}. Since version 2.6.0, it is possible to fit some \dQuote{autocorrelated random-coefficient} models by a syntax consistent with that of random-coefficient terms, \code{<prefix>(<rhs>|.)}. For example, independent Matérn effects can be fitted for males and females by using the syntax \code{Matern(male|.) + Matern(female|.)}, where \code{male} and \code{female} are TRUE/FALSE factors. A numeric variable \code{z} can also be considered, in which case the proper syntax is \code{Matern(0+z|.)}, which represents an autocorrelated random-slope (only) term (or, equivalently, a direct specification of heteroscedasticy of the Matérn random effect). All these effects are achieved by direct control of the elements of the incidence matrix \bold{Z} through the \code{<rhs>} term. By contrast, \code{Matern(z|.)} is not defined. It could mean that a correlation structure between random intercepts and random slopes is to be combined with a Matérn correlation structure, but no way of combining them is yet defined and implemented.  
% cf fns .calc_Z_model_matrix() and .calc_raw_ZA().
% tests: test-spaMM.R and in tests_private/test-syntax-Matern-RHS.R

Since version 2.7.0, the syntax \code{(z-1|.)}, for \bold{numeric} \code{z} only, can also be used to fit some heteroscedastic \emph{non-Gaussian} random effects. For example, a Gamma random-effect term \code{(wei-1|block)} specifies an heteroscedastic Gamma random effect \eqn{u} with constant mean 1 and variance \code{wei^2} \eqn{\lambda}, where \eqn{\lambda} is still the estimated variance parameter. See Details of \code{\link{negbin}} for a possible application. Here, this effect is not implemented through direct control of \bold{Z} (multiplying the elements of an incidence matrix \bold{Z} by \code{wei}), as this would have a different effect on the distribution of the random effect term. \code{(z|.)} is not defined. It could mean that a correlation structure between random intercepts and random slopes for Gamma-distributed random effects is considered, but such correlation structures are not well-specified by their correlation matrix. 

The double-vertical syntax, \code{( rhs || lhs )}, is interpreted as in \code{lme4}. Any such term is immediately converted to \code{( ( 1 | lhs ) + ( 0 + lhs | rhs ) )}, and should be counted as two random effects for all purposes (e.g., for fixing the variances of the random effects). 

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
spaMM is designed to be used through the high-level functions 
\code{\link{fitme}} (the most general fitting function), \code{\link{corrHLfit}}, \code{\link{HLCor}}, \code{\link{HLfit}}, and \code{\link{fixedLRT}}

The \code{test} directory of the package provides many additional examples of spaMM usage beyond those from the formal documentation. 
}

\examples{
## Fit a Poisson GLMM with adjacency (CAR) correlation model
# see ?adjacency for how to fit efficiently such model models 
data("scotlip") ## loads 'scotlip' data frame, but also 'Nmatrix'
HLCor(cases~I(prop.ag/10) +adjacency(1|gridcode)+offset(log(expec)),
          adjMatrix=Nmatrix,family=poisson(),data=scotlip) 


if (spaMM.getOption("example_maxtime")>2.1) {          
 ## Adding a Gamma random effect to fit a negative-binomial response:
 HLCor(cases~I(prop.ag/10) +(1|gridcode)+adjacency(1|gridcode)
                +offset(log(expec)),
          data=scotlip,family=poisson(),rand.family=list(Gamma(log),gaussian()),
          adjMatrix=Nmatrix)
}

 % tested in test-spaMM.R
\dontrun{
 ## fit non-spatial crossed random effects with distinct families
 data("salamander")
 HLfit(cbind(Mate,1-Mate)~1+(1|Female)+(1|Male),family=binomial(),
        rand.family=list(gaussian(),Beta(logit)),data=salamander,HLmethod="ML")
}

 ## Nested effects
 % tested in test-spaMM.R
\dontrun{
 # lmer syntax allowing several degrees of nesting
 HLfit(cbind(Mate,1-Mate)~1+(1|Female/Male),
       family=binomial(),rand.family=Beta(logit),data=salamander,HLmethod="ML")
 # [ also allowed is cbind(Mate,1-Mate)~1+(1|Female)+(1|Male \%in\% Female) ]
}

## fit a non-spatial, Gamma GLMM:
data("wafers")
HLfit(y ~X1*X3+X2*X3+I(X2^2)+(1|batch),family=Gamma(log),
          data=wafers)

## Same with fixed-effects predictor for residual variance 
## ( = structured-dispersion model):
HLfit(y ~X1*X3+X2*X3+I(X2^2)+(1|batch),family=Gamma(log),
          resid.model = ~ X3+I(X3^2) ,data=wafers)

## Random-slope model (mind the output!)        
if (spaMM.getOption("example_maxtime")>1) {          
  HLfit(y~X1+(X2|batch),data=wafers)
}

## fit a GLM (not mixed) with structured dispersion:
HLfit( y ~X1+X2+X1*X3+X2*X3+I(X2^2),family=Gamma(log),
      resid.model = ~ X3+I(X3^2) ,data=wafers)

## Fit of binary data using PQL/L. See ?arabidopsis
\dontrun{
data("arabidopsis")
HLCor(cbind(pos1046738,1-pos1046738)~seasonal+Matern(1|LAT+LONG),
                   ranPars=list(rho=0.129,lambda=4.28,nu=0.291),
                   family=binomial(),HLmethod="PQL/L",data=arabidopsis)
}
}
% ** add an example for fixedLRT ??
\keyword{ package }
\keyword{spatial}
% should add \keyword{smooth} later... and htest ? models ? see RShowDoc("KEYWORDS")
