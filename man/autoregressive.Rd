\name{adjacency}
\alias{adjacency}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fitting autoregressive models
}
\description{
This describes the autoregressive models (AR) considered by \code{spaMM} and how to fit them.
This documentation is mainly a stub, as more implementation of some these models (in particular, of simultaneous AR, i.e., SAR) is to be expected in the future. The only model currently implemented is the \code{adjacency} model (a conditional AR, i.e., CAR). 
}
\details{
The most widely applicable fitting strategy in \code{spaMM} is to fit possibly complex spatial models by numerical maximization over correlation parameters. However, more efficient algorithms are available for CAR and SAR models, and have been widely used in particular in the econometric literature (e.g., LeSage and Pace 2009). Currently such an algorithm is implemented only for one CAR model, the \code{adjacency} model. The taxonomy of AR models is not well settled, but \code{adjacency} should retain its current meaning as a specific CAR model in future versions of \code{spaMM}. 

In the adjacency model, the covariance matrix of random effects \bold{u} can be described as \eqn{\lambda}(\bold{I}\eqn{-\rho} \bold{W}\eqn{)^{-1}} where \bold{W} is the (symmetric) adjacency matrix. If you use \code{corrHLfit}, the default fitting method for the \code{adjacency} model is numerical maximization of the likelihood (or restricted likelihood) as function of the correlation parameter \eqn{\rho}. However, a call to \code{HLCor} fits more efficiently the variance and correlation parameters of random effects (see examples). The ML fits by both methods should be practically equivalent. The REML fits should slightly differ from each other, due to the fact that the REML approximation for GLMMs does not maximize a single likelihood function.
Finally, a call to \code{corrHLfit} with the additional argument 
\code{init.HLfit=list(rho=0)} should be equivalent in speed and result to the \code{HLCor} call.

The efficient method uses the spectral decomposition \bold{W=VDV'} where \bold{D} is a diagonal matrix of eigenvalues \eqn{d_i}, the covariance of \bold{V'u} is 
\eqn{\lambda}(\bold{I}\eqn{-\rho} \bold{D}\eqn{)^{-1}}, which is a diagonal matrix with elements 
\eqn{\lambda_i}=\eqn{\lambda}/(1\eqn{-\rho d_i}). Hence \eqn{1/\lambda_i} is in the linear predictor form \eqn{\alpha}+\eqn{\beta d_i} This can be used to fit \eqn{\lambda} and \eqn{\rho} efficiently. If the efficient fitting method is used, the results are reported as the coefficients \eqn{\alpha} (\code{(Intercept)}) and \eqn{\beta} (\code{adjd}) of the predictor for \eqn{1/\lambda_i}, in addition to the resulting values of \eqn{\rho} and of the common \eqn{\lambda} factor.      
}
\references{
LeSage, J., Pace, R.K. (2009) Introduction to Spatial Econometrics. Chapman & Hall/CRC.
}
\examples{
data(scotlip)
# CAR by Laplace with 'outer' estimation of rho
if (spaMM.getOption("example_maxtime")>10) {          
  corrHLfit(cases~I(prop.ag/10) +adjacency(1|gridcode)+offset(log(scotlip$expec)),
          adjMatrix=Nmatrix,family=poisson(),data=scotlip,HLmethod="ML") ## ~ 10 s.
}

# CAR by Laplace with 'inner' estimation of rho
HLCor(cases~I(prop.ag/10) +adjacency(1|gridcode)+offset(log(scotlip$expec)),
          adjMatrix=Nmatrix,family=poisson(),data=scotlip,HLmethod="ML") ## ~ 1 s.
}