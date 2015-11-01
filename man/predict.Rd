\name{predict}
\alias{predict.HLfit}
\alias{predict}

\title{
  Prediction from a model fit.
}
\description{
  Prediction of the response variable by its expected value obtained as (the inverse link transformation of) the linear predictor (\eqn{\eta}) and more generally for terms of the form \bold{X}[n]\eqn{\beta}+\bold{Z}[n]\bold{v}, for possibly new design matrices \bold{X}[n] and \bold{Z}[n].   
}
\usage{
\method{predict}{HLfit}(object,newdata = newX, newX=NULL, re.form= NULL,
                 variances=list(fixef=FALSE, linPred=FALSE,dispVar=FALSE, 
                                resid=FALSE, sum=FALSE, cov=FALSE),
                 predVar=variances$linPred,residVar=variances$resid,
                 binding = FALSE,...)
}
\arguments{
  \item{object}{
  The return object of an HLfit or similar function.
}
  \item{newdata}{
  \bold{Either} a matrix or data frame containing all required variables for evaluating fixed and random effects, including an offset. If \code{NULL}, the original data are reused.  
  
  \bold{or} a numeric vector, which names (if any) are ignored. This makes it easier to use \code{predict} as an objective function for an 
  optimization procedure such as \code{optim}, which calls the objective function on unnamed vectors. However, one must make sure that the order of elements in the vector is  the order of first occurrence of the variables in the model formula. This order can be checked in the error message returned when calling \code{predict} on a \code{newX} vector of clearly wrong size, e.g. \code{predict(<object>,newdata=numeric(0))}.
}
 \item{newX}{equivalent to newdata, available for back-compatibility} 
 \item{re.form}{
  formula for random effects to include.  By default, it is NULL, in which case all random effects are included. If it is NA, no random effect is included. If it is a formula, only the random effects it contains are retained. The other variance components are removed from both point prediction and \code{variances} calculations. If you want to retain only the spatial effects in the point prediction, but all variances, either use re.form  and add missing variances (on linear predictor scale) manually, or ignore this argument and see Details and Examples for different ways of controlling variances.     
 }
 \item{variances}{
   A list which elements control the computation of different estimated variances. \code{fixef=TRUE} will provide the variances of \bold{X}\eqn{\beta}; \code{linPred=TRUE} will provide the variance of the linear predictor \eqn{\eta} for given dispersion parameters (see Details); \code{resid=TRUE} will provide the residual variances (for Gaussian or Gamma responses). These different variances are returned as attributes \code{"fixefVar"}, \code{"predVar"}, and \code{"residVar"}, respectively. 
   
   If \code{dispVar=TRUE}, \code{"predVar"} may include a component of prediction variance for uncertainty in dispersion parameter estimates. Currently, this option only handles a scalar residual variance (\eqn{\phi}) and a single random effect with a scalar variance (\eqn{\lambda}). 
   
   \code{sum=TRUE} is equivalent to \code{resid=TRUE,linPred=TRUE}: in both cases \code{predVar}, \code{residVar}, and their sum (as attribute \code{"sumVar"}), are returned. If \code{cov=TRUE}, the full covariance matrices are returned for any of the requested terms (except for \code{"residVar"}, as the covariance matrix of the residuals is diagonal).     
 }
 \item{predVar}{
 (for back-compatibility: \code{variances} should now be used)
 \code{predVar=TRUE} corresponds to \code{variances=list(linPred=TRUE)}, and 
 \code{predVar="Cov"} corresponds to \code{variances=list(linPred=TRUE,cov=TRUE)}. 
 }
 \item{residVar}{
 (for back-compatibility: \code{variances} should now be used)
  \code{residVar=TRUE} corresponds to \code{variances=list(resid=TRUE)}.
 }
 \item{binding}{
  If \code{binding} is a character string, the predicted values are bound with the \code{newdata} and the result is returned as a data frame. The predicted values column name is the given \code{binding}, or a name based on it, if the \code{newdata} already include a variable with this name). 
  If \code{binding} is \code{FALSE}, The predicted values are returned as a matrix and the data frame used for prediction is returned as an attribute
  (unless it was \code{NULL}).
 }
\item{\dots}{
  further arguments passed to or from other methods. 
}
}
\value{
  A matrix or data frame (according to the \code{binding} argument), with optionally one or more prediction variance vector or (co)variance matrices as attributes. The further attribute \code{fittedName} contains the binding name, if any.  
}
\details{
  If  \code{newdata} is NULL, \code{predict} returns the fitted responses, including random effects, from the object. 
  Otherwise it computes new predictions including random effects as far as possible.   
  For spatial random effects it constructs a correlation matrix \bold{C} between new locations and locations in the original fit. Then it infers the random effects in the new locations as   \bold{C} (\bold{L}'\eqn{)^{-1}} \bold{v} (see \code{\link{spaMM}} for notation).   For non-spatial random effects, it checks whether any group (i.e., level of a random effect) in the new data was represented in the original data, and it adds the inferred random effect for this group to the prediction for individuals in this group. 

  \code{fixefVar} is the (co)variance of \bold{X}\eqn{\beta} (or \bold{X}[n]\eqn{\beta}), deduced from the asymptotic covariance matrix of \eqn{\beta} estimates. 

  \code{predVar} is the prediction (co)variance of \eqn{\eta}=bold{X}\eqn{\beta}+\bold{ZL}\bold{v} (see \code{\link{HLfit}} Details for notation), or more generally of \bold{X}[n]\eqn{\beta}+\bold{Z}[n]\bold{L}[n]\bold{v}, by default computed for given dispersion parameters. It assumes that the covariance matrix of \eqn{\beta} and \bold{v} estimates is the inverse of the expected Hessian matrix (for given dispersion parameters) of the augmented linear model for \eqn{\beta} and \bold{v}. It thus takes into account the joint uncertainty in estimation of \eqn{\beta} and prediction of \bold{v}.

For \bold{prediction covariance} with a new \bold{Z}[n], it matters whether a single or multiple new levels are used: see Examples.      

If \code{variances$dispVar} is \code{TRUE}, it may also include a term accounting for uncertainty in \eqn{\phi} and \eqn{\lambda} (Booth and Hobert 1998, eq. 19). 

  For models with non-Gaussian response, the prediction covariance of the response is approximated by the prediction covariance of the linear predictor, 
  pre- and post-multiplied by \eqn{\partial\mu/\partial\eta}. 
  
  These variance calculations are approximate except for LMMs, and ignore uncertainties in spatial correlation parameters, 
  and thus cannot be garanteed to give accurate results.  

In the \bold{point prediction} of the linear predictor, the unconditional expected value of \eqn{u} is assigned to the realizations of \eqn{u} for unobserved levels of non-spatial random effects (it is zero in GLMMs but not for non-gaussian random effects), and the inferred value of \eqn{u} is assigned in all other cases. Corresponding values of \eqn{v} are then deduced. This computation yields the classical \dQuote{BLUP} or empirical Bayes predictor in LMMs, but otherwise it may yield less well characterized predictors, where \dQuote{unconditional} \eqn{v} may not be its expected value when the \code{rand.family} link is not identity. 
}
\references{
  Booth, J.G., Hobert, J.P. (1998) Standard errors of prediction in generalized linear mixed models. J. Am. Stat. Assoc. 93: 262-272. 
}

\examples{
data(blackcap)
fitobject <- corrHLfit(migStatus ~ 1 + Matern(1|latitude+longitude),data=blackcap,
                       ranFix=list(nu=4,rho=0.4,phi=0.05))
predict(fitobject)
predict(fitobject,blackcap) ## same computation, different format 
getDistMat(fitobject)

#### multiple controls of prediction variances
## (1) fit with an additional random effect
grouped <- cbind(blackcap,grp=c(rep(1,7),rep(2,7))) 
fitobject <- corrHLfit(migStatus ~ 1 +  (1|grp) +Matern(1|latitude+longitude),
                       data=grouped,  ranFix=list(nu=4,rho=0.4,phi=0.05))

## (2) re.form usage to remove a random effect from point prediction and variances: 
predict(fitobject,re.form= ~ 1 +  Matern(1|latitude+longitude))

## (3) comparison of covariance matrices for two types of new data
moregroups <- grouped[1:5,]
rownames(moregroups) <- paste("newloc",1:5,sep="")
moregroups$grp <- rep(3,5) ## all new data belong to an unobserved third group 
cov1 <- attr(predict(fitobject,newdata=moregroups,
                     variances=list(linPred=TRUE,cov=TRUE)),"predVar")
moregroups$grp <- 3:7 ## all new data belong to distinct unobserved groups
cov2 <- attr(predict(fitobject,newdata=moregroups,
                     variances=list(linPred=TRUE,cov=TRUE)),"predVar")
cov1-cov2 ## the expected off-diagonal covariance due to the common group in the first fit.

\dontrun{
## Effects of numerically singular correlation matrix C:
fitobject <- corrHLfit(migStatus ~ 1 + Matern(1|latitude+longitude),data=blackcap,
                       ranFix=list(nu=10,rho=0.001)) ## numerically singular C
predict(fitobject) ## predicted mu computed as X beta + L v 
predict(fitobject,newdata=blackcap) ## predicted mu computed as X beta + C %*% solve(t(L),v) 

## prediction with new X and Z
if(require("rsae", quietly = TRUE)) {
  data(landsat)
  fitobject <- HLfit(HACorn ~ PixelsCorn + PixelsSoybeans + (1|CountyName),
                     data=landsat[-33,],HLmethod="ML")
  newXandZ <- unique(data.frame(PixelsCorn=landsat$MeanPixelsCorn,
                                PixelsSoybeans=landsat$MeanPixelsSoybeans,
                                CountyName=landsat$CountyName))
  predict(fitobject,newdata=newXandZ,variances = list(linPred=TRUE,dispVar=TRUE))
}

}
}
%- \keyword{ print }
