\name{numInfo}
\alias{numInfo}
\title{
Information matrix
}
\description{
Computes the information matrix for (ideally) all parameters for main response model, using numerical derivatives, that is the matrix of second derivatives of negative log likelihood. Residual dispersion parameters are not handled beyond scalar \code{phi} values.  % ___F I X M E___
The default value of the \code{which} argument shows all classes of parameters that should be handled, including random-effect parameters (\code{lambda}, \code{\link{ranCoefs}}, \code{\link{corrPars}}, and \code{\link{hyper}}), residual dispersion parameters (\code{phi}, \code{NB_shape} for \code{\link{negbin1}} and \code{\link{negbin2}}, and \code{beta_prec} for \code{\link{beta_resp}}), and fixed-effect coefficients (\code{beta}).
}
\usage{
numInfo(fitobject, transf = FALSE, 
        which = c("lambda", "ranCoefs", "corrPars", "hyper", 
                  "phi", "NB_shape", "beta_prec", "beta"), 
        verbose=FALSE, ...)
}
\arguments{
  \item{fitobject}{
  Fit object returned by a \pkg{spaMM} fitting function.
}
  \item{transf}{
  Whether to perform internal computations on a transformed scale (but computation on 
  transformed scale may be implemented for fewer classes of models than default computation). 
}
  \item{which}{
  character vector: set of parameters with respect to which derivatives are to be computed. 
}
  \item{verbose}{
Boolean: whether to print (as a list) the estimates of the parameters for which the Hessian will be computed, additional information about possibly ignored parameters, possible misure of REML fits, and a (sort of) progress bar if the procedure is expected to last more than a few seconds.
}
  \item{\dots}{Arguments passed to \code{numDeriv::hessian} and \code{numDeriv::grad}.}
}
\value{
A matrix
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
