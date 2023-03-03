\name{ARp}
\alias{ARp}
\alias{ARMA}
\title{
Random effect with AR(p) (autoregressive of order p) or ARMA(p,q) structure.
}
\description{
These times-series correlation models can be declared as correlation models for random effect. 
The AR(p) model is here parametrized by the \bold{partial} correlation coefficients of the levels of the random effect, \{U_t\}, corr(U_s,U_t|U_(s+1),...,U_(t-1)), with valid values in the hypercube ]-1,,1[^p (Barndorff-Nielsen and Schou, 1973). 
In the autoregressive-moving average ARMA(p,q) model, the AR part is parametrized in the same way. AR parameters are named \code{"p1"}, \code{"p2"}..., and MA parameters are named \code{"q1"}, \code{"q2"}... .

Implementation of the AR(p) model uses the sparsity of the inverse covariance matrix. In the ARMA(p,q) model, neither the covariance nor its inverse are sparse, so fits are expected to be more time- and memory-consuming. 


}
\usage{
# corrFamily constructors:
ARp(p=1L, fixed=NULL, corr=TRUE, tpar=1/(1+seq(p)))
ARMA(p=1L, q=1L, fixed=NULL, tpar=c(1/(1+seq_len(p)),1/(1+seq_len(q))))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{Integer: order of the autoregressive process.}
  \item{q}{Integer: order of the moving-average process.}
  \item{tpar}{
  Numeric vector: template values of the \bold{partial coefficient coefficients} of the autoregressive process, and the traditional coefficients of the moving-average processe, in this order. The \code{tpar} vector must always have full length, even when some parameters are fixed.}
  \item{fixed}{NULL or numeric vector, to fix the parameters of this model.}
  \item{corr}{For development purposes, better ignored in normal use.}

}
\references{
Barndorff-Nielsen 0. and Schou G., 1973 On the parametrization of autoregressive models by partial autocorrelations. J. Multivariate Analysis 3: 408-419. \doi{10.1016/0047-259X(73)90030-4}
}
\value{
The \code{ARp} and \code{ARMA} functions return a \code{\link{corrFamily}} descriptor, hence a \code{list} including element \code{Cf}, a function returning, for given ARMA or AR parameters, the correlation matrix for \code{ARMA}, or its \bold{inverse} for \code{ARp}.

The fitted correlation matrix can be extracted from a fit object, as for any autocorrelated random effect, by \code{Corr(}<fit object>\code{)[[}<random-effect index>\code{]]}. 
}
\examples{
if (spaMM.getOption("example_maxtime")>2) {
ts <- data.frame(lh=lh,time=seq(48)) ## using 'lh' data from 'stats' package

## Default 'tpar' => AR1 model 
#
(ARpfit <-  fitme(lh ~ 1 + ARp(1|time), data=ts, method="REML"))
#
## which is equivalent to
#
(AR1fit <- fitme(lh ~ 1 +AR1(1|time), data=ts, method="REML"))

## AR(3) model 
#
(AR3fit <- fitme(lh ~ 1 + ARp(1|time, p=3), data=ts, method="REML"))

## Same but with fixed 2-lag partial autocorrelation 
#
(AR3fix <- fitme(lh ~ 1 + ARp(1|time, p=3, fixed=c(p2=0)), data=ts, method="REML"))
#      
# The fit should be statistically equivalent to
#
(AR3_fix <- fitme(lh ~ 1 + ARp(1|time, p=3), data=ts, method="REML",
                 fixed=list(corrPars=list("1"=c(p2=0)))))
#                 
# with subtle differences in the structure of the fit objects:
#
get_ranPars(AR3fix)$corrPars     # p2 was not a parameter of the model
get_ranPars(AR3_fix)$corrPars    # p2 was a fixed parameter of the model
#
# get_fittefPars() expectedly ignores 'p2' whichever way it was fixed.    


## Same as 'AR3fix' but with an additional MA(1) component
#
(ARMAfit <- fitme(lh ~ 1 + ARMA(1|time, p=3, q=1, fixed=c(p2=0)), 
                  data=ts, method="REML"))
}


}
\keyword{ ts }
