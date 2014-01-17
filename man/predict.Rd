\name{predict}
%- \alias{predict.HL}
\alias{predict}
\alias{predict.HLfit} 
\alias{predict.HLCor}
\alias{predict.corrHLfit}

\title{
  Prediction from a model fit.
}
\description{
  Predictions of the response variable in different locations and/or different values of the predictor variables for fixed effects, based on random effects inferred by HLCor / corrHLfit.  
}
\usage{
%- \method{predict}{HL}(object,newX = NULL,coeffs=NULL,...)
\method{predict}{HLfit}(object,newX = NULL,...)
\method{predict}{HLCor}(object,newX = NULL,coeffs=NULL,...)
\method{predict}{corrHLfit}(object,newX = NULL,coeffs=NULL,...)
}
\arguments{
  \item{object}{
  The return object of an HLfit / HLCor / corrHLfit call.
}
  \item{newX}{
  A data frame containing new locations and values of the predictor variables for fixed effects and of any variable required to compute an offset. Its column names must match those of the originally fitted data. If \code{NULL}, the original data are reused. If not \code{NULL}, then the original fit should not be based on an adjacency matrix. 
}
 \item{coeffs}{
  Precomputed coefficients for the prediction (see Details).
 }
\item{\dots}{
  further arguments passed to or from other methods. 
}
}
\value{
  A data frame including the locations and the expected value (\code{fv=}\eqn{\mu}) of the response variable. 
}
\details{
  If  \code{newX} is NULL, \code{predict} only returns the fitted responses from the object.
  Otherwise it constructs a correlation matrix \bold{C} between new locations and locations in the original fit. Then it infers the random effects in the new locations as 
  \bold{C} [t(\bold{L})]^{-1} \bold{v} (see \code{\link{spaMM}} for notation). The fixed effects and offset (if any) are added. If the predictor is used many times, it may be useful to precompute [t(\bold{L})]^{-1} \bold{v} and to provide this vector through the \code{coeffs} argument (see Examples).  

}

\examples{
data(blackcap)
fitobject <- corrHLfit(migStatus ~ 1 + Matern(1|latitude+longitude),data=blackcap,
                       ranFix=list(nu=4,rho=0.4,phi=0.05))
predict(fitobject)
## same result using precomputed 'coeffs':
coeffs <- with(fitobject$hlfit,solve(t(attr(predictor,"LMatrix")),v_h))
predict(fitobject,blackcap,coeffs=coeffs)

\dontrun{
## Effects of numerically singular correlation matrix C:
fitobject <- corrHLfit(migStatus ~ 1 + Matern(1|latitude+longitude),data=blackcap,
                       ranFix=list(nu=10,rho=0.001)) ## numerically singular C
predict(fitobject) ## predicted mu computed as X beta + L v 
predict(fitobject,newX=blackcap) ## predicted mu computed as X beta + C.solve(t(L),v) 
}
}
%- \keyword{ print }
