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
Diverse autoregressive (AR) models are implemented in spaMM. This documentation describe the \code{adjacency} model (a conditional AR, i.e., CAR), and the \code{AR1} model for time series. Other documentation deals with more or less distantly related models: \code{\link{ARp}} for more general AR(p) and ARMA(p,q) models for time series, and \code{\link{IMRF}} and \code{\link{MaternIMRFa}} for mesh-based approximations of geostatistical models.

An AR1 random effect is specified as \code{AR1(1|<grouping factor>)}. It describes correlations between realizations of the random effect for (typically) successive time-steps by a correlation \eqn{\phi}, denoted \code{ARphi} in function calls. Nested AR1 effects can be specified by a nested grouping factor, as in \code{AR1(1|<time index> \%in\% <nesting factor>)}.

A CAR random effect is specified as \code{adjacency(1|<grouping factor>)}. The correlations among levels of the random effect form a 
    matrix (\bold{I}\eqn{-\rho} \code{adjMatrix}\eqn{)^{-1}}, in terms of an \code{adjMatrix} matrix which must be provided, and of the scalar \eqn{\rho}, denoted \code{rho} in function calls. 
    The rows and columns of \code{adjMatrix} must have names matching those of levels of the random effect \bold{or else} be ordered as increasing values of the levels of the geographic location index specifying the spatial random effect. For example, if the model formula 
    is\cr 
    \code{y ~ adjacency(1|geo.loc)} and \code{<data>$geo.loc} is 2,4,3,1,... the first row/column of the matrix refers to \code{geo.loc}=1, i.e. to the fourth row of the data. 
}
\details{
Efficient algorithms for CAR models have been widely discussed in particular in the econometric literature (e.g., LeSage and Pace 2009), but these models are not necessarily recommended for irregular lattices (see Wall, 2004 and Martellosio, 2012 for some insights on the implications of autoregressive models). 

In \bold{CAR} models, the covariance matrix of random effects \bold{u} can be described as \eqn{\lambda}(\bold{I}\eqn{-\rho} \bold{W}\eqn{)^{-1}} where \bold{W} is the (symmetric) adjacency matrix. \code{HLCor} uses the spectral decomposition of the adjacency matrix, written as \bold{W=VDV'} where \bold{D} is a diagonal matrix of eigenvalues \eqn{d_i}. The covariance of \bold{V'u} is 
\eqn{\lambda}(\bold{I}\eqn{-\rho} \bold{D}\eqn{)^{-1}}, which is a diagonal matrix with elements 
\eqn{\lambda_i}=\eqn{\lambda}/(1\eqn{-\rho d_i}). Hence \eqn{1/\lambda_i} is in the linear predictor form \eqn{\alpha}+\eqn{\beta d_i} This can be used to fit \eqn{\lambda} and \eqn{\rho} efficiently. A call to \code{corrHLfit} with the additional argument 
\code{init.HLfit=list(rho=0)} should be equivalent in speed and result to the \code{HLCor} call. 

This is fast for small datasets (as in the example below) but more generic maximization algorithms may be preferable for large ones. It is suggested to use \code{fitme} generally unless one has a large number of small data sets to analyze. A call to \code{fitme} or \code{corrHLfit} without that initial value does not use the spectral decomposition. It performs numerical maximization of the likelihood (or restricted likelihood) as function of the correlation parameter \eqn{\rho}. The choice of fitting function may slightly impact the results. The ML fits by \code{corrHLfit} and \code{HLCor} should be practically equivalent. The REML fits should slightly differ from each other, due to the fact that the REML approximation for GLMMs does not maximize a single likelihood function. 

If \code{HLCor} is used, the results are reported as the coefficients \eqn{\alpha} (\code{(Intercept)}) and \eqn{\beta} (\code{adjd}) of the predictor for \eqn{1/\lambda_i}, in addition to the resulting values of \eqn{\rho} and of the common \eqn{\lambda} factor.    

Different fits may also differ in using or not algorithms that exploit the sparsity of the precision matrix of the autoregressive random effect. By default, spaMM tends to select sparse-precision algorithms for large datasets and large (i.e. many-level) random effects (details are complex).
However, for \bold{AR1} models, the dimension of the implied precision matrix is determined by the extreme values of grouping factor (typically interpreted as a time index), as all intermediate values must be considered. Then, the correlation-based algorithms may be more efficient if only a few levels are present in the data, as only a small correlation matrix is required in that case.    
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
fitme(lh ~ 1 +AR1(1|time), data=ts, method="REML")
# With fixed parameters:
# HLCor(lh ~ 1 +AR1(1|time), data=ts, ranPars=list(ARphi=0.5,lambda=0.25,phi=0.001))

##### CAR random effect:
data("scotlip")
# CAR by Laplace with 'outer' estimation of rho
if (spaMM.getOption("example_maxtime")>0.8) {          
  fitme(cases ~ I(prop.ag/10)+adjacency(1|gridcode)+offset(log(expec)),
          adjMatrix=Nmatrix, family=poisson(), data=scotlip) 
}

# CAR by Laplace with 'inner' estimation of rho
HLCor(cases ~ I(prop.ag/10)+adjacency(1|gridcode)+offset(log(expec)),
          adjMatrix=Nmatrix, family=poisson(), data=scotlip, method="ML")
}
\keyword{ models }% use one of  RShowDoc("KEYWORDS")
\keyword{ spatial }% __ONLY ONE__ keyword per line
\keyword{ ts }% __ONLY ONE__ keyword per line
