\name{numInfo}
\alias{numInfo}
\alias{print.singeigs}
\title{
Information matrix
}
\description{
Computes by numerical derivation the observed information matrix for (ideally) all parameters for mean response model, that is, the matrix of second derivatives of negative log likelihood. 
The default value of the \code{which} argument shows all classes of parameters that should be handled, including random-effect parameters (\code{lambda}, \code{\link{ranCoefs}}, \code{\link{corrPars}}, and \code{\link{hyper}}), residual dispersion parameters (\code{phi}, \code{NB_shape} for \code{\link{negbin1}} and \code{\link{negbin2}}, and \code{beta_prec} for \code{\link{beta_resp}} and \code{\link{betabin}}), and fixed-effect coefficients (\code{beta}).

Model fits including a \code{\link{phi-resid.model}} are not fully handled, in two ways: the information matrix does not include their parameters; and if the residual dispersion model include random effects, there is good reason for the \code{numInfo} calculation to detect that the fit has not maximized marginal likelihood with respect to most parameters.

}
\usage{
numInfo(fitobject, transf = FALSE, which = NULL, check_deriv = TRUE,
        sing=1e-05, verbose=FALSE, refit_hacks=list(), ...)
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
Boolean: whether to print (as a list) the estimates of the parameters for which the Hessian will be computed, additional information about possibly ignored parameters, possible misuse of REML fits, and a (sort of) progress bar if the procedure is expected to last more than a few seconds.
}
  \item{refit_hacks}{
  list of arguments; its name anticipates that it might allow hazardous manipulations in a later version of \pkg{spaMM}. But currently only the innocuous element \code{verbose} of the list will be taken into account. Notably, \code{refit_hacks=list(verbose=c(TRACE=TRUE))} can be used to give information on parameter values used in the computation of numerical derivatives.
  }
  \item{\dots}{Arguments passed to \code{numDeriv::hessian} and \code{numDeriv::grad}.}
}
\details{
The computation of a second derivatives is not necessarily meaningful if a first derivative does not vanish at the parameter estimates. This may occur in particular when the objective function (say, marginal likelihood) is maximized at a boundary of the parameter space (say, at zero for \code{lambda} estimates). Further, for such values at a boundary, only one-sided derivatives can be computed, and this is not handled by \code{numDeriv::hessian}. So, some checks may be requested to detect non-zero gradients and parameters estimated at their boundaries. The boundary checks are currently performed for \code{lambda} and \code{ranCoefs} estimates, if \code{check_deriv} is set to TRUE or to NULL. Other parameters are not (yet) checked, so \code{numInfo} may sometimes fails when such other parameter estimates are at a boundary. If \code{check_deriv} is set to TRUE, an additional systematic check of the gradient with respect to all target parameters is performed. 
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
(SEs <- sqrt(diag(solve(numinfo))))
#
# => beta SEs here equal to conditional SEs shown by fit summary.
# Other SEs can be compared to the approximate ones 
# for log(phi) and log(lambda), given by 
#
# update(lmmfit, control=list(refit=TRUE))
#
# => 1118*0.5289 and 10840*0.1024

data("blackcap")
maternfit <- fitme(migStatus ~ means+ Matern(1|longitude+latitude),data=blackcap) 
numInfo(maternfit)
}
