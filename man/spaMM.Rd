\name{spaMM}
% ----- \Rdversion{1.1}
\alias{spaMM}

\docType{package}
\encoding{latin1}

\title{Inference in spatial GLMMs (and a few other models)}

\description{spaMM is used to fit spatial GLMMs.
The package can fit a wider class of non-spatial models (see Examples) but has not been optimized for that purpose. 
It fits models which linear predictor is of the form \eqn{y ~ }linear predictor for fixed effects    +   spatially correlated random effect    +    \code{\link{offset}}.

The random effects are either Gaussian (which defines GLMMs), or other distributions (which defines the wider class of hierarchical GLMs).  
}

\details{
   The models are fit via Laplace approximations to compute (1) the marginal likelihood with respect to random effects and (2) the restricted likelihoods (as in REML), 
     i.e conditional likelihoods with respect to fixed effect estimates. 

  The package has been thoroughly tested only for analysis of spatial GLMMs, 
  where the random effects are Gaussian. The alternative of Beta-distributed, Gamma and inverse-Gamma distributed random effects 
  have been implemented as detailed in the \code{\link{HLfit}} documentation. 
  They have not been thoroughly tested (though see again \code{\link{HLfit}} for a valid example of the latter). 
  They are not recommended for spatial models.    

  The variance of random effects is denoted \eqn{\lambda} and that of residual error is denoted \eqn{\phi}.
  A fixed-effects linear predictor for \eqn{\phi} can be considered (see Examples). Fixed effects are described in the standard form \bold{X}\eqn{\beta}
  where \bold{X} is the design matrix of fixed effects and \eqn{\beta} a vector of fixed effect parameters. The correlated random effects can be described 
  as \bold{ZLv}, where \bold{ZL} is the product of the \bold{Z} and \bold{L} matrices described below and forms the \dQuote{design matrix} for random effects,
  and \bold{v} is a vector of independent and identically distributed random effects.
  The package makes it easy to incorporate spatial correlations in the \bold{ZL} matrix, following either a 
  Matern correlation model or an adjacency matrix model (see \code{\link{HLCor}}). 

  The structure of the model fitted can be described by the following steps. first, independent and identically distributed (iid) random effects \eqn{u} are drawn from some distribution;
  second, a transformation \bold{v}\eqn{=f}(\bold{u}) is applied (\bold{v} elements are still iid); third, correlated random effects are obtained as \bold{Lv} 
  where \bold{L} is the \dQuote{square root} of a correlation matrix (this may be meaningful only for Gaussian random effects). 
  Finally, a matrix \bold{Z} allows to specify how the correlated random effects
  affect the response values, so that the design matrix for \bold{v} is the product \bold{ZL}. 
 
}
\author{\enc{Fran?ois}{Francois} Rousset and Jean-Baptiste Ferdy}
\references{
Lee, Y., Nelder, J. A. and Pawitan, Y. (2006). Generalised linear models with random effects: unified analysis via
h-likelihood. Chapman & Hall: London.
}
\seealso{
This is designed to be used through the documented high-level functions 
\code{\link{corrHLfit}}, \code{\link{HLCor}}, \code{\link{HLfit}},  \code{\link{Predictor}}, \code{\link{fixedLRT}}
}

\examples{
## Fit a Poisson GLMM with adjacency correlation model

\dontrun{
data(scotlip)
corrHLfit(Predictor(cases~I(prop.ag/10) +adjacency(1|gridcode),offset=log(scotlip$expec)),
          data=scotlip,family=poisson(),
          adjMatrix=Nmatrix,lower=list(rho=0),upper=list(rho=0.1745))
}

## fit a non-spatial, Gamma-logNormal model
data(wafers)
HLfit(predictor= y ~X1*X3+X2*X3+I(X2^2)+(1|batch),family=Gamma(log),
          resid.predictor = "phi" ~ 1,data=wafers)

## Same with fixed-effects predictor for residual variance.
HLfit(predictor= y ~X1*X3+X2*X3+I(X2^2)+(1|batch),family=Gamma(log),
          resid.predictor = "phi" ~ X3+I(X3^2) ,data=wafers)

## Fit of binary data using PQL/L. See ?arabidopsis
\dontrun{
data(arabidopsis)
HLCor(cbind(pos1046738,1-pos1046738)~seasonal+Matern(1|LAT+LONG),
                   ranPars=list(rho=0.129,lambda=4.28,nu=0.291),
                   family=binomial(),HLmethod="PQL/L",data=arabidopsis)
}
}
% ** add an example for fixedLRT ??
\keyword{ package }
\keyword{spatial}
% should add \keyword{smooth} later... and htest ? models ? see RShowDoc("KEYWORDS")
