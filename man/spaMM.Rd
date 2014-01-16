\name{spaMM}
% ----- \Rdversion{1.1}
\alias{spaMM}

\docType{package}

\title{Inference in spatial GLMMs (and other non-spatial mixed models)}

\description{spaMM has been designed to fit spatial mixed models.
This package can fit a wider class of non-spatial mixed models (see Examples), but has not been optimized for that purpose. 
The random effects are either Gaussian (which defines GLMMs), or other distributions (which defines the wider class of hierarchical GLMs), or simply absent (which makes a GLM).  
}

\details{
   The models are fit via Laplace approximations to compute (1) the marginal likelihood with respect to random effects and (2) the restricted likelihoods (as in REML), 
     i.e conditional likelihoods with respect to fixed effect estimates. 

  The package has been thoroughly tested mainly for analysis of spatial GLMMs, 
  where the random effects are Gaussian. The alternative of Beta-distributed, Gamma and inverse-Gamma distributed random effects 
  have been implemented as detailed in the \code{\link{HLfit}} documentation and checked against literature results, 
  but they are not recommended for spatial models.    

  The variance of random effects is denoted \eqn{\lambda} and that of residual error is denoted \eqn{\phi}.
  A fixed-effects linear predictor for \eqn{\phi}, modeling heteroskedasticity, can be considered (see Examples). Fixed effects are described in the standard form \bold{X}\eqn{\beta}
  where \bold{X} is the design matrix of fixed effects and \eqn{\beta} a vector of fixed effect parameters. 
  
  The high-level functions of the package make it easy to incorporate spatial correlations in the model, following either a 
  Matérn correlation model or an adjacency matrix model. 
  
  For more complex applications, the structure of the random effects can be described by the following steps. 
  First, independent and  identically distributed (iid) random effects \eqn{u} are drawn from some distribution;
  second, a transformation \bold{v}\eqn{=f}(\bold{u}) is applied (\bold{v} elements are still iid); third, correlated random effects are obtained as \bold{Mv}, where  
  \bold{M} can describe spatial correlation between observed locations, block effects (or repeated observations in given locations), and correlations involving unobserved locations. See Details in \code{\link{Predictor}} for further details.   
 
}
\author{\enc{François}{Francois} Rousset and Jean-Baptiste Ferdy}
\references{
Lee, Y., Nelder, J. A. and Pawitan, Y. (2006). Generalised linear models with random effects: unified analysis via
h-likelihood. Chapman & Hall: London.
}
\seealso{
spaMM is designed to be used through the high-level functions 
\code{\link{corrHLfit}}, \code{\link{HLCor}}, \code{\link{HLfit}},  \code{\link{Predictor}}, \code{\link{fixedLRT}}
}

\examples{
## Fit a Poisson GLMM with adjacency correlation model

\dontrun{
data(scotlip) ## loads scotlip data frame, and Nmatrix
corrHLfit(cases~I(prop.ag/10) +adjacency(1|gridcode)+offset(log(scotlip$expec)),
          data=scotlip,family=poisson(),
          adjMatrix=Nmatrix,lower=list(rho=0),upper=list(rho=0.1745))
}

## fit a non-spatial, Beta-binomial model:
data(Loaloa)
idx  <-as.factor(seq(nrow(Loaloa)))
fLoaloa <- cbind(Loaloa,idx)
HLfit(cbind(npos,ntot-npos)~elev1+elev2+elev3+elev4+maxNDVI1+seNDVI
                   +(1|idx),HLmethod="HL(1,1)",
            data=fLoaloa,family=binomial(),rand.family=Beta()) 

## fit a non-spatial, Gamma GLMM:
data(wafers)
HLfit(y ~X1*X3+X2*X3+I(X2^2)+(1|batch),family=Gamma(log),
          data=wafers)

## Same with fixed-effects predictor for residual variance 
## ( = structured-dispersion model):
HLfit(y ~X1*X3+X2*X3+I(X2^2)+(1|batch),family=Gamma(log),
          resid.formula = ~ X3+I(X3^2) ,data=wafers)

## fit a GLM (not mixed) with structured dispersion:
HLfit( y ~X1+X2+X1*X3+X2*X3+I(X2^2),family=Gamma(log),
      resid.formula = ~ X3+I(X3^2) ,data=wafers)

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
