\name{autoregressive}
\alias{autoregressive}
\alias{adjacency}
\alias{CAR}
\alias{AR1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fitting autoregressive models
}
\description{
Two autoregressive(AR) models are currently implemented: the \code{adjacency} model (a conditional AR, i.e., CAR), and the \code{AR1} model for time series. Implementation of further models (in particular, of simultaneous AR, i.e., SAR) is to be expected in the future.
Efficient algorithms for CAR models have been widely discussed in particular in the econometric literature (e.g., LeSage and Pace 2009), but these models are not necessarily recommended for irregular lattices (see Wall, 2004 and Martellosio, 2012 for some insights on the implications of autoregressive models). The fastest method in spaMM for large data sets is implemented in the \code{fitme} function. For small data sets (as in the example below), \code{HLCor} may be fastest. It is suggested to use \code{fitme} generally unless one has a large number of small data sets to analyze. For non-LMMs, \code{corrHLfit(*,HLmethod="PQL/L")} can be quite fast (but does not return an ML or REML fit).   

An AR1 random effect is specified as \code{AR1(1|<grouping factor>)}. It describes correlations between realizations of the random effect for (typically) successive time-steps by a correlation \eqn{\phi}, denoted \code{ARphi} in function calls. Nested AR1 effects can be specified by a nested grouping factor, as in \code{AR1(1|<time index> \%in\% <nesting factor>)}.

A CAR random effect is specified as \code{adjacency(1|<grouping factor>)}. The correlations among levels of the random effect form a 
    matrix (\bold{I}\eqn{-\rho} \code{adjMatrix}\eqn{)^{-1}}, in terms of an \code{adjMatrix} matrix which must be provided, and of the scalar \eqn{\rho}, denoted \code{rho} in function calls. 
    The rows and columns of \code{adjMatrix} must be ordered as increasing values of the levels of the geographic 
    location index specifying the spatial random effect. For example, if the model formula 
    is\cr 
    \code{y ~ adjacency(1|geo.loc)} and \code{<data>$geo.loc} is 2,4,3,1,... the first row/column of the matrix refers to \code{geo.loc}=1, i.e. to the fourth row of the data. 
}
\details{
%The taxonomy of AR models is not well settled, but \code{adjacency} should retain its current meaning as a specific CAR model in future versions of \code{spaMM}.

For \bold{AR1} models, and large data sest, the fitting functions by default select methods that exploits the sparsity of the precision matrix of the random effects. the dimension of the implied precision matrix is determined by the extreme values of grouping factor (typically interpreted as a time index), as all intermediate values must be considered. Thus, the precision matrix may be quite large even if few levels are represented in the data.   

For \bold{CAR} models, different fitting strateges may be used:

A call to \code{HLCor} uses the spectral decomposition of the adjacency matrix as further detailed below. This is fast for small datasets but fitme may be preferable otherwise.

A call to \code{corrHLfit} with the additional argument 
\code{init.HLfit=list(rho=0)} should be equivalent in speed and result to the \code{HLCor} call.

A call to \code{corrHLfit} without this argument does not use the spectral decomposition. It performs a generic numerical maximization of the likelihood (or restricted likelihood) as function of the correlation parameter \eqn{\rho}. The ML fits by \code{corrHLfit} and \code{HLCor} should be practically equivalent. The REML fits should slightly differ from each other, due to the fact that the REML approximation for GLMMs does not maximize a single likelihood function. 

In the adjacency model, the covariance matrix of random effects \bold{u} can be described as \eqn{\lambda}(\bold{I}\eqn{-\rho} \bold{W}\eqn{)^{-1}} where \bold{W} is the (symmetric) adjacency matrix. \code{HLCor} uses the spectral decomposition of the adjacency matrix, written as bold{W=VDV'} where \bold{D} is a diagonal matrix of eigenvalues \eqn{d_i}. The covariance of \bold{V'u} is 
\eqn{\lambda}(\bold{I}\eqn{-\rho} \bold{D}\eqn{)^{-1}}, which is a diagonal matrix with elements 
\eqn{\lambda_i}=\eqn{\lambda}/(1\eqn{-\rho d_i}). Hence \eqn{1/\lambda_i} is in the linear predictor form \eqn{\alpha}+\eqn{\beta d_i} This can be used to fit \eqn{\lambda} and \eqn{\rho} efficiently. If \code{HLCor} is used, the results are reported as the coefficients \eqn{\alpha} (\code{(Intercept)}) and \eqn{\beta} (\code{adjd}) of the predictor for \eqn{1/\lambda_i}, in addition to the resulting values of \eqn{\rho} and of the common \eqn{\lambda} factor.      
}
\references{
LeSage, J., Pace, R.K. (2009) Introduction to Spatial Econometrics. Chapman & Hall/CRC.

Martellosio, F. (2012) The correlation structure of spatial autoregressions, Econometric Theory 28, 1373-1391.

Wall M.M. (2004) A close look at the spatial structure implied by
the CAR and SAR models: Journal of Statistical Planning and Inference 121: 311-324.
}
\examples{
##### AR1 random effect:
ts <- data.frame(lh=lh,time=seq(48)) ## using 'lh' data from stats package
HLCor(lh~1 +AR1(1|time), data=ts, ranPars=list(ARphi=0.5,lambda=0.25,phi=0.001))


##### CAR random effect:
data("scotlip")
# CAR by Laplace with 'outer' estimation of rho
if (spaMM.getOption("example_maxtime")>0.7) {          
  corrHLfit(cases~I(prop.ag/10) +adjacency(1|gridcode)+offset(log(expec)),
          adjMatrix=Nmatrix,family=poisson(),data=scotlip,HLmethod="ML") 
}
if (spaMM.getOption("example_maxtime")>0.8) {          
  fitme(cases~I(prop.ag/10) +adjacency(1|gridcode)+offset(log(expec)),
          adjMatrix=Nmatrix,family=poisson(),data=scotlip) 
}

# CAR by Laplace with 'inner' estimation of rho
HLCor(cases~I(prop.ag/10) +adjacency(1|gridcode)+offset(log(expec)),
          adjMatrix=Nmatrix,family=poisson(),data=scotlip,HLmethod="ML")
}