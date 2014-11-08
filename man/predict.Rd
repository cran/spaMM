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
                 object,newX = NULL,coeffs=NULL,re.form= NULL,
                 variances=list(fixef=FALSE,ranef=FALSE,resid=FALSE,sum=FALSE,cov=FALSE),
                 predVar=variances$ranef,residVar=variances$resid,
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
 \item{re.form}{
  formula for random effects to include.  If NULL, include all
          random effects; if NA, include no random effects
 }
 \item{variances}{
   A list which elements control whether to compute different estimated variances (and their sum). \code{fixef=TRUE} will provide the variances of \bold{X}\eqn{\beta}; \code{ranef=TRUE} will provide the prediction variance of the random effects; \code{resid=TRUE} will provide the residual variances (for Gaussian or Gamma responses). These different variances are returned as attributes \code{"fixefVar"}, \code{"predVar"} and \code{"residVar"}, respectively. If \code{sum=TRUE}, all components are computed. If more than one of the components is computed, their sum is returned as attribute \code{"sumVar"}. If \code{cov=TRUE}, the full covariance matrices are returned for any of the requested terms (except for \code{"residVar"}, as the covariance matrix of the residuals is diagonal). See Details for how \code{"predVar"} is computed.    
 }
 \item{predVar}{
 (for back-compatibility: \code{variances} should now be used)
 \code{predVar=TRUE} corresponds to \code{variances=list(ranef=TRUE)}, and 
 \code{predVar="Cov"} corresponds to \code{variances=list(ranef=TRUE,cov=TRUE)}. 
 }
 \item{residVar}{
 (for back-compatibility: \code{variances} should now be used)
  \code{residVar=TRUE} corresponds to \code{variances=list(resid=TRUE)}.
 }
 \item{binding}{
  If \code{binding} is a valid variable name for a data frame, the predicted values are bound (under the given name) 
  with the data frame used for prediction and the resulting frame is returned.
  If \code{binding} is \code{FALSE}, The predicted values are returned as a matrix and the data frame used for prediction is returned as an attribute
  (unless it was \code{NULL}).
 }
\item{\dots}{
  further arguments passed to or from other methods. 
}
}
\value{
  A matrix or data frame (according to the \code{binding} argument), with optionally one or more prediction variance vector or (co)variance matrices as attributes.  
}
\details{
  If  \code{newX} is NULL, \code{predict} only returns the fitted responses, including random effects, from the object. 
  Otherwise it computes new predictions including random effects as far as possible.   
  For spatial random effects it constructs a correlation matrix \bold{C} between new locations and locations in the original fit. Then it infers the random effects in the new locations as   \bold{C} (\bold{L}'\eqn{)^{-1}} \bold{v} (see \code{\link{spaMM}} for notation). If the predictor is used many times, it may be useful to precompute (\bold{L}'\eqn{)^{-1}} \bold{v} and to provide this vector through the \code{coeffs} argument (see Examples).  
  For non-spatial random effects, it checks whether any group (i.e., level of a random effect) in the new data was represented in the original data, and it adds the inferred random effect for this group to the prediction for individuals in this group.
  
  The \bold{prediction variance} is the variance of the linear predictor (\eqn{\eta}). \code{"predVar"} is the prediction variance of the random effect terms in (\eqn{\eta}). It takes into account the uncertainty in estimation of \eqn{\beta}, and is computed as described in Gotway and Wolfinger (2003) based on earlier works for LMMs. 
  
  Unobserved levels of non-spatial random effects are handled as follows. In the \bold{point prediction} of the linear predictor, 
  the expected value of \eqn{u} is assigned to the realizations of \eqn{u} for unobserved groups (this value is 0 in LMMs). Corresponding realizations of \eqn{v} are then deduced using the link function(s) for the random effects (the identity link in LMMs). The same computation is performed in all other models, for good or bad. For \bold{prediction covariance}, it matters whether a single or multiple new levels are used: see Examples.    
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
predict(fitobject,coeffs=coeffs,variances=list(sum=TRUE)) -> pf
attr(pf,"sumVar")


###### handling of unobserved groups
## (1) fit with an additional random effect
grouped <- cbind(blackcap,grp=c(rep(1,7),rep(2,7))) 
fitobject <- corrHLfit(migStatus ~ 1 +  (1|grp) +Matern(1|latitude+longitude),
                       data=grouped,  ranFix=list(nu=4,rho=0.4,phi=0.05))
## (2) comparison of covariance matrices for two types of new data
moregroups <- grouped[1:5,]
rownames(moregroups) <- paste("newloc",1:5,sep="")
moregroups$grp <- rep(3,5) ## all new data belong to an unobserved third group 
cov1 <- attr(predict(fitobject,newX=moregroups,
                     variances=list(ranef=TRUE,cov=TRUE)),"predVar")
moregroups$grp <- 3:7 ## all new data belong to distinct unobserved groups
cov2 <- attr(predict(fitobject,newX=moregroups,
                     variances=list(ranef=TRUE,cov=TRUE)),"predVar")
cov1-cov2 ## the expected off-diagonal covariance due to the common group in the first fit.

\dontrun{
## Effects of numerically singular correlation matrix C:
fitobject <- corrHLfit(migStatus ~ 1 + Matern(1|latitude+longitude),data=blackcap,
                       ranFix=list(nu=10,rho=0.001)) ## numerically singular C
predict(fitobject) ## predicted mu computed as X beta + L v 
predict(fitobject,newX=blackcap) ## predicted mu computed as X beta + C %*% solve(t(L),v) 
}
}
%- \keyword{ print }
