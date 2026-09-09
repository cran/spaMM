\name{numInfo}
\alias{numInfo}
\alias{print.singeigs}
\title{
Information matrix
}
\description{
Computes by numerical derivation the observed information matrix for (ideally) all parameters for mean response model, that is, the matrix of second derivatives of negative log likelihood. 
The default value of the \code{which} argument shows all classes of parameters that should be handled (except for PQL/L fits: see Details), including random-effect parameters (\code{lambda}, \code{\link{ranCoefs}}, \code{\link{corrPars}}, and \code{\link{hyper}}), residual dispersion parameters (\code{phi}, \code{NB_shape} for \code{\link{negbin1}} and \code{\link{negbin2}}, and \code{beta_prec} for \code{\link{beta_resp}} and \code{\link{betabin}}), and fixed-effect coefficients (\code{beta}).

Approximate covariance matrices of parameters estimates can be deduced from the inverse of the information matrix. Inverses of blocks of this matrix provide asymptotic covariance matrices given the values of (nuisance) parameters excluded from a given block: these covariances calculations ignore the additional variance stemming from the uncertainty in the nuisance parameter values. For example, the inverse of the block for fixed-effect coefficients should contain variances corresponding to the standard errors traditionally reported in tables of fixed effects (see Examples). 

The function calls algorithms from \pkg{numDeriv} and share their limitations. Notably, they may request a re-fit of the model with some parameters fixed at values that are out of the allowed ranges. This can be controlled by the \code{method.args} argument, passed through the \dots to these algorithms (see Examples). 

PQL/L fits raise specific problems, as do other fits by methods maximizing two different likelihood approximations for different subsets of parameters (see Details).

Model fits including a \code{\link{phi-resid.model}} are not fully handled, in two ways: the information matrix does not include their parameters; and if the residual dispersion model include random effects, there is good reason for the \code{numInfo} calculation to detect that the fit has not maximized marginal likelihood with respect to most parameters.



}
\usage{
numInfo(fitobject, transf = FALSE, which = NULL, check_deriv = TRUE,
        sing=1e-05, verbose=FALSE, refit_hacks=list(), attrs=NULL, 
        return.="", ...)
}
\arguments{
  \item{fitobject}{
  Fit object returned by a \pkg{spaMM} fitting function.
}
  \item{transf}{
  Whether to perform internal computations on a transformed scale (but computation on 
  transformed scale may be implemented for fewer classes of models than default computation). 
}
  \item{which}{ NULL, or character vector giving the sets of parameters with respect to which derivatives are to be computed.
  The NULL default is equivalent to \code{c("lambda", "ranCoefs", "corrPars", "hyper", "phi", "NB_shape", "beta_prec", "beta")} for ML fits,
  and to the same except \code{"beta"} (fixed effects) for REML fits.
}
  \item{check_deriv}{
  Boolean; whether to perform some checks for possible problems (see Details). 
  }
  \item{sing}{
  numeric value, or \code{FALSE}; if it is a nonzero numeric value, eigenvalues of the matrix are checked and values smaller than \code{sing}
  are highlighted in output (see Value). This will highlight nearly-singular information matrices, but also those with large negative eigenvalues.
  }
  \item{verbose}{
Boolean: whether to print (as a list) the estimates of the parameters for which the information matrix will be computed, additional information about possibly ignored parameters, possible misuse of REML fits, and a (sort of) progress bar if the procedure is expected to last more than a few seconds.
}
  \item{refit_hacks}{
  list of arguments; its name anticipates that it might allow hazardous manipulations in a later version of \pkg{spaMM}. But currently only the innocuous element \code{verbose} of the list will be taken into account. Notably, \code{refit_hacks=list(verbose=c(TRACE=TRUE))} can be used to give information on parameter values used in the computation of numerical derivatives.
  }
  \item{attrs}{NULL or vector of character strings, for names of optional attributes to be added to the return value; in particular, \code{"df"} for a data frame of the parameter points and their likelihood values used to compute the information matrix.}
  \item{return.}{character string. If set to \code{"grad"}, the function returns the numerical gradient of the objective function instead of the numerical information matrix.}
  \item{\dots}{Arguments passed to \code{\link[numDeriv]{hessian}} and \code{\link[numDeriv]{grad}}, or functions with a similar interface. In particular, \code{method.args} can be used as shown in the Examples.}
}
\details{
The computation of a second derivatives is not necessarily meaningful if a first derivative does not vanish at the parameter estimates. This may occur in particular when the objective function (say, marginal likelihood) is maximized at a boundary of the parameter space (say, at zero for \code{lambda} estimates). Further, for such values at a boundary, only one-sided derivatives can be computed, and this is not handled by \code{numDeriv::hessian}. So, some checks may be requested to detect non-zero gradients and parameters estimated at their boundaries. The boundary checks are currently performed for \code{lambda} and \code{ranCoefs} estimates, if \code{check_deriv} is set to TRUE or to NULL. Other parameters are not (yet) checked, so \code{numInfo} may sometimes fails when such other parameter estimates are at a boundary. If \code{check_deriv} is set to TRUE, an additional systematic check of the gradient with respect to all target parameters is performed. 

For PQL/L fits, the gradient of the log-likelihood approximation used to infer random-effect parameters (\eqn{\theta}, say) is not zero at the fixed-effect estimates, since the fixed-effect coefficients \eqn{\beta} are then estimated by maximizing the distinct h-likelihood. Further, trying to compute the full information matrix (including fixed-effect coefficients) in this case means that fixed effects are fixed to they estimates \eqn{\hat{\beta}} when \eqn{\theta} is varied, while \eqn{\hat{\theta}} maximize the likelihood approximation only when \eqn{\beta} is refitted for any given \eqn{\theta}. Thus, the gradient for \eqn{\theta} with \eqn{\beta} fixed is also not zero. For these reasons, the only information matrix that can be safely computed excludes the fixed effects, as in the latter case the \eqn{\beta} are refitted for any given \eqn{\theta}.

}
\value{
NULL or a matrix.

NULL is returned if no parameter is found with respect to which a numerical information \dQuote{should} be computed (where what should be done depends on the \code{which} and \code{check_derivs} arguments). 

Otherwise a matrix is returned, with an \code{eigvals} attribute if \code{sing} was non-zero. This attribute is a numeric vector of eigenvalues of the matrix. If some eigenvalue(s) were lower than \code{sing}, the vector has additional class \code{"singeigs"} so that its printing is controlled by an ad-hoc \code{print.singeigs} method highlighting the small eigenvalues. 
}
\examples{
data("wafers")
lmmfit <- fitme(y ~X1+X2+X1*X3+X2*X3+I(X2^2)+(1|batch),data=wafers)
numinfo <- numInfo(lmmfit)
covmat <- solve(numinfo[-(1:2),-(1:2)])
(cond.SEs <- sqrt(diag(covmat))) # = SEs for fixed effects in the summary()

(SEs <- sqrt(diag(solve(numinfo)))) 
# => SEs[-(1:2)] here, *but not generally*, equivalent to cond.SEs.
# SEs[1:2] are those of lambda and phi; they can be compared 
# to the more approximate ones, for log(phi) and log(lambda), given by 
# > update(lmmfit, control=list(refit=TRUE))
# These SEs translate to 1118*0.5289 ~ 591.3... for lambda and 
# 10840*0.1024 ~ 1110 for phi, while numInfo() gives 735.9 and 1120, respectively.

data("blackcap")
maternfit <- fitme(migStatus ~ means+ Matern(1|longitude+latitude),data=blackcap) 
numInfo(maternfit)

\dontrun{
# Using 'method.args':
set.seed(123)
fx <- as.vector(na.omit(stats::filter(rnorm(400L), rep(0.025,40))))
tsl <- length(fx)
ts1 <- data.frame(y=fx+rnorm(tsl, sd=0.1),time=seq(tsl))
ts2 <- data.frame(y=fx+rnorm(tsl, sd=0.1),time=seq(tsl))
toydata <- rbind(ts1,ts2)
(fts <- fitme(y ~ 1 +AR1(1|time), data=toydata))
# numInfo(fts) # Fails as a value of AR1 correlation parameter > 1 is tried 
numInfo(fts, method.args=list(d=2e-5)) # OK
}

}
