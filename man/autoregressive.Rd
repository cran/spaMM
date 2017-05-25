\name{adjacency}
\alias{adjacency}
\alias{CAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fitting autoregressive models
}
\description{
The only autoregressive model (AR) currently implemented is the \code{adjacency} model (a conditional AR, i.e., CAR). Implementation of further models (in particular, of simultaneous AR, i.e., SAR) is to be expected in the future.
Efficient algorithms have been widely discussed in particular in the econometric literature (e.g., LeSage and Pace 2009). The fastest method in spaMM for large data sets is implemented in the \code{fitme} function. For small data sets (as in the example below), \code{HLCor} may be fastest. It is suggested to use \code{fitme} generally unless one has a large number of small data sets to analyze. For LMMs, \code{corrHLfit(*,HLmethod="PQL/L")} can be quite fast (but does not return an ML or REML fit).   
}
\details{
%The taxonomy of AR models is not well settled, but \code{adjacency} should retain its current meaning as a specific CAR model in future versions of \code{spaMM}.
For large data sest, \code{fitme} selects methods that exploits the sparsity of the precision matrix of the random effects. 

A call to \code{HLCor} uses the spectral decomposition of the adjacency matrix as further detailed below. 

A call to \code{corrHLfit} with the additional argument 
\code{init.HLfit=list(rho=0)} should be equivalent in speed and result to the \code{HLCor} call.

A call to \code{corrHLfit} without this argument does not use the spectral decomposition. It performs a generic numerical maximization of the likelihood (or restricted likelihood) as function of the correlation parameter \eqn{\rho}. The ML fits by \code{corrHLfit} and \code{HLCor} should be practically equivalent. The REML fits should slightly differ from each other, due to the fact that the REML approximation for GLMMs does not maximize a single likelihood function. 

In the adjacency model, the covariance matrix of random effects \bold{u} can be described as \eqn{\lambda}(\bold{I}\eqn{-\rho} \bold{W}\eqn{)^{-1}} where \bold{W} is the (symmetric) adjacency matrix. \code{HLCor} uses the spectral decomposition of the adjacency matrix, written as bold{W=VDV'} where \bold{D} is a diagonal matrix of eigenvalues \eqn{d_i}. The covariance of \bold{V'u} is 
\eqn{\lambda}(\bold{I}\eqn{-\rho} \bold{D}\eqn{)^{-1}}, which is a diagonal matrix with elements 
\eqn{\lambda_i}=\eqn{\lambda}/(1\eqn{-\rho d_i}). Hence \eqn{1/\lambda_i} is in the linear predictor form \eqn{\alpha}+\eqn{\beta d_i} This can be used to fit \eqn{\lambda} and \eqn{\rho} efficiently. If \code{HLCor} is used, the results are reported as the coefficients \eqn{\alpha} (\code{(Intercept)}) and \eqn{\beta} (\code{adjd}) of the predictor for \eqn{1/\lambda_i}, in addition to the resulting values of \eqn{\rho} and of the common \eqn{\lambda} factor.      
}
\references{
LeSage, J., Pace, R.K. (2009) Introduction to Spatial Econometrics. Chapman & Hall/CRC.
}
\examples{
data(scotlip)
# CAR by Laplace with 'outer' estimation of rho
if (spaMM.getOption("example_maxtime")>2.2) {          
  corrHLfit(cases~I(prop.ag/10) +adjacency(1|gridcode)+offset(log(expec)),
          adjMatrix=Nmatrix,family=poisson(),data=scotlip,HLmethod="ML") 
}
if (spaMM.getOption("example_maxtime")>1.6) {          
  fitme(cases~I(prop.ag/10) +adjacency(1|gridcode)+offset(log(expec)),
          adjMatrix=Nmatrix,family=poisson(),data=scotlip) 
}

# CAR by Laplace with 'inner' estimation of rho
HLCor(cases~I(prop.ag/10) +adjacency(1|gridcode)+offset(log(expec)),
          adjMatrix=Nmatrix,family=poisson(),data=scotlip,HLmethod="ML")
}