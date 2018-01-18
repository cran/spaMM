\name{spaMM}
% ----- \Rdversion{1.1}
\alias{spaMM}
\alias{spaMM-package}

\docType{package}

\title{Inference in mixed models, in particular spatial GLMMs}

\description{Fits a range of mixed models, including those with spatially correlated random effects. The random effects are either Gaussian (which defines GLMMs), or other distributions (which defines the wider class of hierarchical GLMs), or simply absent (which makes a GLM).  
}

\details{

  The standard response families \code{gaussian}, \code{binomial}, \code{poisson}, and \code{Gamma} are handled, as well as negative binomial (see \code{\link{negbin}}), zero-truncated poisson and negative binomial, and Conway-Maxwell-Poisson response (see \code{\link{Tpoisson}}, \code{\link{Tnegbin}} and \code{\link{COMPoisson}}). A \code{multi} family look-alike is also available for \code{\link{multinomial}} response, with some constraints. 
  
  The package fits models including several nested or crossed random effects, one of which can be following a Matérn correlation model (see \code{\link{Matern}}), or an adjacency matrix model (see \code{\link{adjacency}}), or an \code{AR1} model, or a gven \code{corrMatrix}. GLMMs and HGLMs are fit via Laplace approximations for (1) the marginal likelihood with respect to random effects and (2) the restricted likelihood (as in REML), i.e the likelihood of random effect parameters given the fixed effect estimates. 
  
  The variance(s) of random effects (\eqn{u}) is (are) denoted \eqn{\lambda} (\code{lambda} in input and output). The variance parameter of residual error is denoted \eqn{\phi} (\code{phi}): this is the residual variance for gaussian response, but for Gamma-distributed response, the residual variance is \eqn{\phi}\eqn{\mu^2} where \eqn{\mu} is expected response. A fixed-effects linear predictor for \eqn{\phi}, modeling heteroscedasticity, can be considered (see Examples). Fixed effects are described in the standard form \bold{X}\eqn{\beta} where \bold{X} is the design matrix of fixed effects and \eqn{\beta} (\code{beta}) is a vector of fixed effect parameters. 

  The structure of the random effects can generally be described by the following steps. First, independent and  identically distributed (iid) random effects \eqn{u} are drawn from one of the following distributions: gaussian, Beta-distributed, Gamma and   inverse-Gamma distributed random effects, implemented as detailed in the \code{\link{HLfit}} documentation. Second, a transformation \bold{v}\eqn{=f}(\bold{u}) is applied (\bold{v} elements are still iid). Third, correlated random effects are obtained as \bold{Mv}, where the matrix \bold{M} can describe spatial correlation between observed locations, block effects (or repeated observations in given locations), and correlations involving unobserved locations. See Details in \code{\link{Predictor}} for the general form of \bold{M}. In most cases \bold{M} is determined from the model formula, but it can also be input directly (e.g., to describe genetic correlations).     

  The package has been extensively tested mainly for analysis of spatial GLMMs (Rousset and Ferdy 2014 and ongoing work), where the random effects are Gaussian. Other models have been checked against literature results and a few simulations.
 
}
\author{\enc{François}{Francois} Rousset and Jean-Baptiste Ferdy.

The syntax of formula terms mostly matches the one in the \code{lme4} package, so bits of code for parsing formulas has been recycled in \code{spaMM} from various versions of \code{lme4}. 
}
\references{
Lee, Y., Nelder, J. A. and Pawitan, Y. (2006). Generalized linear models with random effects: unified analysis via
h-likelihood. Chapman & Hall: London.

Rousset F., Ferdy, J.-B. (2014) Testing environmental and genetic effects in the presence of spatial autocorrelation. Ecography, 37: 781-790.
    \url{http://dx.doi.org/10.1111/ecog.00566}

}
\seealso{
spaMM is designed to be used through the high-level functions 
\code{\link{corrHLfit}}, \code{\link{HLCor}}, \code{\link{HLfit}}, and \code{\link{fixedLRT}}
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
