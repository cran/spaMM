\name{predict}
\alias{predict.HLfit}
\alias{predict}

\title{
  Prediction from a model fit.
}
\description{
  Predictions of the response variable, based on given values of the predictor variables for fixed effects, and/or 
  on predicted random effects.   
}
\usage{
\method{predict}{HLfit}(
                 object,newX = NULL,coeffs=NULL,
                 predVar=FALSE, re.form= NULL,
                 binding = if(is.vector(newX)) {FALSE} else {"fitted"},...)
}
\arguments{
  \item{object}{
  The return object of an HLfit or similar function.
}
  \item{newX}{
  \bold{Either} a matrix or data frame containing all required variables for evaluating fixed and random effects, including an offset. If \code{NULL}, the original data are reused.  
  
  \bold{or} a numeric vector, which names (if any) are ignored. This makes it easier to use \code{predict} as an objective function for an 
  optimization procedure such as \code{optim}, which calls the objective function on unnamed vectors. However, one must make sure that the order of elements in the vector is  the order of first occurrence of the variables in the model formula. This order can be checked in the error message returned when calling \code{predict} on a \code{newX} vector of clearly wrong size, e.g. \code{predict(<object>,newX=numeric(0))}.
}
 \item{coeffs}{
  Precomputed coefficients for the prediction (see Details).
 }
 \item{predVar}{
  Either a boolean or \code{"Cov"}. If this is \code{TRUE}, the prediction variances are computed. If this is \code{"Cov"}, the full
  covariance matrix is computed. In both cases, the result is included as attribute \code{"predVar"} of the return value. See Details for how prediction variances are computed. 
 }
 \item{re.form}{
  formula for random effects to include.  If NULL, include all
          random effects; if NA, include no random effects
 }
 \item{binding}{
  If \code{binding} is a valid variable name for a data frame, the predicted values are bound (under the given name) 
  with the date frame used for prediction and the resulting frame is returned.
  If \code{binding} is \code{FALSE}, The predicted values are returned as a matrix and the data frame used for prediction is returned as an attribute
  (unless it was \code{NULL}).
 }
\item{\dots}{
  further arguments passed to or from other methods. 
}
}
\value{
  A matrix or data frame (according to the \code{binding} argument), with optionally a prediction variance vector or (co)variance matrix as attribute.  
}
\details{
  If  \code{newX} is NULL, \code{predict} only returns the fitted responses, including random effects, from the object. 
  Otherwise it computes new predictions including random effects as far as possible.   
  For spatial random effects it constructs a correlation matrix \bold{C} between new locations and locations in the original fit. Then it infers the random effects in the new locations as   \bold{C} (\bold{L}'\eqn{)^{-1}} \bold{v} (see \code{\link{spaMM}} for notation). If the predictor is used many times, it may be useful to precompute (\bold{L}'\eqn{)^{-1}} \bold{v} and to provide this vector through the \code{coeffs} argument (see Examples).  
  For non-spatial random effects, it checks whether any group (i.e., level of a random effect) in the new data was represented in the original data, and it adds the inferred random effect for this group to the prediction for individuals in this group.
  
  The \bold{prediction variance} is the variance of the linear predictor (\eqn{\eta}). It accounts for the uncertainty in prediction of the random effects. In addition, it accounts for uncertainty in estimation of \eqn{\beta}. It is computed as described in Gotway and Wolfinger (2003). 
}
\references{
  Gotway, C.A., Wolfinger, R.D. (2003) Spatial prediction of counts and rates. Statistics in Medicine 22: 1415-1432. 
}

\examples{
data(blackcap)
fitobject <- corrHLfit(migStatus ~ 1 + Matern(1|latitude+longitude),data=blackcap,
                       ranFix=list(nu=4,rho=0.4,phi=0.05))
predict(fitobject)

predict(fitobject,blackcap) ## same computation, different format 

## same result using precomputed 'coeffs':
coeffs <- predictionCoeffs(fitobject) ## using dedicated extractor function
predict(fitobject,coeffs=coeffs,predVar=TRUE)
attr(.Last.value,"predVar")

\dontrun{
## Effects of numerically singular correlation matrix C:
fitobject <- corrHLfit(migStatus ~ 1 + Matern(1|latitude+longitude),data=blackcap,
                       ranFix=list(nu=10,rho=0.001)) ## numerically singular C
predict(fitobject) ## predicted mu computed as X beta + L v 
predict(fitobject,newX=blackcap) ## predicted mu computed as X beta + C %*% solve(t(L),v) 
}
}
%- \keyword{ print }