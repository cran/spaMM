\name{predict}
\alias{predict.HLfit}
\alias{predict}
\alias{get_fixefVar}
\alias{get_predVar}
\alias{get_residVar}
\alias{get_respVar}
\alias{get_intervals}
\alias{intervals}
\alias{get_predCov_var_fix}
\alias{preprocess_fix_corr}

\title{
  Prediction from a model fit.
}
\description{
  Prediction of the response variable by its expected value obtained as (the inverse link transformation of) the linear predictor (\eqn{\eta}) and more generally for terms of the form \bold{X}_n\eqn{\beta}+\bold{Z}_n\bold{L}\bold{v}, for new design matrices \bold{X}_n and \bold{Z}_n. Various components of prediction variances and predictions intervals can also be computed using \code{predict}. The \code{get_}... functions are convenient extractors for such components. \code{get_predCov_var_fix} extracts a block of a prediction covariance matrix. It was conceived for the specific purpose of computing the spatial prediction covariances between two ``new'' sets of geographic locations, without computing the full covariance matrix for both the new locations and the original (fitted) locations. \code{preprocess_fix_corr} was conceived to perform this covariance computation more efficiently when one of the two sets of new locations is fixed while the other varies.  
}
\usage{
\method{predict}{HLfit}(object, newdata = newX, newX = NULL, re.form = NULL,
                          variances=list(), binding = FALSE, intervals = NULL,
                          level = 0.95, blockSize = 2000L, ...)
get_predCov_var_fix(object, newdata = NULL, fix_X_ZAC.object, fixdata, 
                   variances=list(disp=TRUE,residVar=FALSE), ...)    
preprocess_fix_corr(object, fixdata, re.form = NULL,
                   variances=list(residVar=FALSE))                                
get_fixefVar(...)
get_predVar(...)
get_residVar(...)
get_respVar(...)
get_intervals(...)
}
\arguments{
  \item{object}{
  The return object of fitting functions \code{HLfit,corrHLfit,HLCor}... returning an object inheriting from \code{HLfit} class.
}
  \item{newdata}{
  \bold{Either} NULL, a matrix or data frame, or a numeric vector. 
  
  If \code{NULL}, the original data are reused. Otherwise, all variables required to evaluate model formulas must be included. Which variables are required may depend on other arguments: see \dQuote{prediction with given phi's} example, also illustrating the syntax when formulas include an offset.  
  
  \bold{or} a numeric vector, which names (if any) are ignored. This makes it easier to use \code{predict} as an objective function for an 
  optimization procedure such as \code{optim}, which calls the objective function on unnamed vectors. However, one must make sure that the order of elements in the vector is  the order of first occurrence of the variables in the model formula. This order can be checked in the error message returned when calling \code{predict} on a \code{newX} vector of clearly wrong size, e.g. \code{predict(<object>,newdata=numeric(0))}.
}
 \item{newX}{equivalent to newdata, available for back-compatibility} 
 \item{re.form}{
  formula for random effects to include.  By default, it is NULL, in which case all random effects are included. If it is NA, no random effect is included. If it is a formula, only the random effects it contains are retained. The other variance components are removed from both point prediction and \code{variances} calculations. If you want to retain only the spatial effects in the point prediction, but all variances, either use re.form  and add missing variances (on linear predictor scale) manually, or ignore this argument and see Details and Examples for different ways of controlling variances.     
 }
 \item{variances}{
 

   A list which elements control the computation of different estimated variances. In particular,
   \code{list(linPred=TRUE,disp=TRUE)} is suitable for uncertainty in point prediction.   

   \code{predict} can return four components of prediction variance: 
   \code{fixefVar}, \code{predVar}, \code{residVar} and \code{respVar}, detailed below. 
   They are all returned as attributes of the point predictions. By default, 
   each component is a vector of variances. However, if \code{variances$cov=TRUE}, 
   a covariance matrix is returned when applicable (i.e. not for \code{"residVar"}).
   
   \code{fixefVar} is the (co)variance of fixed effects (\bold{X}\eqn{\beta}) 
   due to uncertainty in \eqn{\beta}. It is called by \code{variances$fixefVar=TRUE}.
   
   \code{predVar} is the (co)variance of the linear predictor \eqn{\eta}. 
   It is the most common measure of uncertainty in point prediction.
   It accounts for uncertainty in fixed effects (\bold{X}\eqn{\beta}) and 
   random effects (\bold{ZLv}) for given dispersion parameters (see Details), 
   but it can also accounts for uncertainty in dispersion parameters 
   (\eqn{\lambda} and \eqn{phi}) estimates if \code{variances$disp=TRUE}, 
   for models in which the effect of uncertainty in dispersion parameters can be computed.  
   Currently, this effect can be computed for a scalar residual variance (\eqn{\phi}) 
   and a single random effect with a scalar variance (\eqn{\lambda}). 
   \code{variances$predVar=TRUE} will return the sum of the two components, if available; 
   otherwise it returns only the (co)variance for given \eqn{\lambda} and \eqn{\phi}. 
   The latter component can be requested by \code{variances$linPred=TRUE}. 
   
   \code{residVar} provides the residual variances (for Gaussian or Gamma responses). 
   It is called by \code{variances$residVar=TRUE}.
  
  \code{respVar} is the variance of the response (see Details).
  It is called by \code{variances$respVar=TRUE}.
    
  Calling for one (co)variance implies that some of its components may be also returned.    

 }
 \item{intervals}{
   NULL or character string or vector of strings. Provides prediction intervals with nominal level \code{level}, deduced from the given prediction variance term, e.g. \code{intervals="predVar"}. Currently only intervals from \code{fixefVar} and \code{predVar} (and for LMMs \code{respVar} including the residual variance) may have a probabilistic meaning. Intervals returned in other cases are (currently) meaningless. 
 }
 \item{level}{Coverage of the intervals.}
 \item{binding}{
  If \code{binding} is a character string, the predicted values are bound with the \code{newdata} and the result is returned as a data frame. The predicted values column name is the given \code{binding}, or a name based on it if the \code{newdata} already include a variable with this name. 
  If \code{binding} is \code{FALSE}, The predicted values are returned as a one-column matrix and the data frame used for prediction is returned as an attribute
  (unless it was \code{NULL}). If \code{binding} is \code{NA}, a vector is returned, without the previous attributes.
 }
 \item{fixdata}{A data frame describing reference data which covariances with variable \code{newdata} may be requested.}
 \item{fix_X_ZAC.object}{The return value of calling \code{preprocess_fix_corr} (see trivial Example). This is a more efficient way of providing information about the \code{fixdata} for repeated calls to \code{get_predCov_var_fix} with variable \code{newdata}.}
 \item{blockSize}{
   Mainly for development purposes. For original or new data with many rows, it may be more efficient to split these data in small blocks, and this gives the maximum number or rows of the blocks. However, this will be ignored if a prediction covariance matrix is requested.    
 }
\item{\dots}{
  further arguments passed to or from other methods. For the \code{get_}... functions, they are passed to \code{predict}. 
}
}
\value{
  For \code{predict}, a matrix or data frame (according to the \code{binding} argument), with optional attributes \code{frame}, \code{intervals}, \code{predVar}, \code{fixefVar}, \code{residVar}, and/or \code{respVar}, the last four holding one or more variance vector or covariance matrices. The further attribute \code{fittedName} contains the binding name, if any.  
  
  The \code{get_}... extractor functions call \code{predict} and extract from its result the attribute implied by the name of the extractor. By default, \code{get_intervals} will return prediction intervals using \code{predVar}.
}
\details{
  If  \code{newdata} is NULL, \code{predict} returns the fitted responses, including random effects, from the object. 
  Otherwise it computes new predictions including random effects as far as possible.   
  For spatial random effects it constructs a correlation matrix \bold{C} between new locations and locations in the original fit. Then it infers the random effects in the new locations as   \bold{C} (\bold{L}'\eqn{)^{-1}} \bold{v} (see \code{\link{spaMM}} for notation).   For non-spatial random effects, it checks whether any group (i.e., level of a random effect) in the new data was represented in the original data, and it adds the inferred random effect for this group to the prediction for individuals in this group. 

  \code{fixefVar} is the (co)variance of \bold{X}\eqn{\beta} (or \bold{X}_n\eqn{\beta}), deduced from the asymptotic covariance matrix of \eqn{\beta} estimates. 

  \code{predVar} is the prediction (co)variance of \eqn{\eta}=\bold{X}\eqn{\beta}+\bold{Z}\bold{v} (see \code{\link{HLfit}} Details for notation), or more generally of \bold{X}_n\eqn{\beta}+\bold{Z}_n\bold{L}\bold{v}, by default computed for given dispersion parameters. 
  
  For levels of the random effects present in the original data, \code{predVar} computation 
%assumes that the covariance matrix of \eqn{\beta} and \bold{v} estimates is the inverse of the expected Hessian matrix (for given dispersion parameters) of the augmented linear model for \eqn{\beta} and \bold{v}. It thus 
takes into account the joint uncertainty in estimation of \eqn{\beta} and prediction of \bold{v}. 
% No mention of augmented model in doc !
  
For new levels of the random effects, \code{predVar} computation additionally takes into account uncertainty in prediction of \bold{v} for these new levels. For \bold{prediction covariance} with a new \bold{Z}_n, it matters whether a single or multiple new levels are used: see Examples.      

If \code{variances$disp} is \code{TRUE}, prediction variance may also include a term accounting for uncertainty in \eqn{\phi} and \eqn{\lambda}, computed following Booth and Hobert (1998, eq. 19). This computation is currently implemented for models with a single random effect, and ignore uncertainties in spatial correlation parameters.

\code{respVar} is the sum of \code{predVar} (pre- and post-multiplied by \eqn{\partial\mu/\partial\eta}{d\mu/d\eta} for models with non-identity link) and of \code{residVar}. 
  
These variance calculations are approximate except for LMMs, and cannot be guaranteed to give accurate results.  

In the \bold{point prediction} of the linear predictor, the unconditional expected value of \eqn{u} is assigned to the realizations of \eqn{u} for unobserved levels of non-spatial random effects (it is zero in GLMMs but not for non-gaussian random effects), and the inferred value of \eqn{u} is assigned in all other cases. Corresponding values of \eqn{v} are then deduced. This computation yields the classical \dQuote{BLUP} or empirical Bayes predictor in LMMs, but otherwise it may yield less well characterized predictors, where \dQuote{unconditional} \eqn{v} may not be its expected value when the \code{rand.family} link is not identity. 

\bold{Intervals} computations use the relevant variance estimates plugged in a Gaussian approximation, except for the simple linear model where it uses Student's \emph{t} distribution.   

}
\references{
  Booth, J.G., Hobert, J.P. (1998) Standard errors of prediction in generalized linear mixed models. J. Am. Stat. Assoc. 93: 262-272. 
}

\examples{
data("blackcap")
fitobject <- corrHLfit(migStatus ~ 1 + Matern(1|latitude+longitude),data=blackcap,
                       ranFix=list(nu=4,rho=0.4,phi=0.05))
predict(fitobject)
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
cov1 <- get_predVar(fitobject,newdata=moregroups,
                     variances=list(linPred=TRUE,cov=TRUE))
moregroups$grp <- 3:7 ## all new data belong to distinct unobserved groups
cov2 <- get_predVar(fitobject,newdata=moregroups,
                     variances=list(linPred=TRUE,cov=TRUE))
cov1-cov2 ## the expected off-diagonal covariance due to the common group in the first fit.

\dontrun{
## prediction with distinct given phi's in different locations:
varphi <- cbind(blackcap,logphi=runif(14))
vphifit <- corrHLfit(migStatus ~ 1 + Matern(1|latitude+longitude), 
                     resid.model = list(formula=~0+offset(logphi)),
                     data=varphi,  ranFix=list(nu=4,rho=0.4))
# for respVar computation, one needs the resid.model formula to specify phi:
get_respVar(vphifit,newdata=data.frame(latitude=1,longitude=1,logphi=1))
# for predVar computation, phi is not needed 
#     (and could have been specified through ranFix):  
get_predVar(vphifit,newdata=data.frame(latitude=1,longitude=1))                     

## Effects of numerically singular correlation matrix C:
fitobject <- corrHLfit(migStatus ~ 1 + Matern(1|latitude+longitude),data=blackcap,
                       ranFix=list(nu=10,rho=0.001)) ## numerically singular C
predict(fitobject) ## predicted mu computed as X beta + L v 
predict(fitobject,newdata=blackcap) ## predicted mu computed as X beta + C %*% solve(t(L),v) 
#
fix_X_ZAC.object <- preprocess_fix_corr(fitobject,fixdata=blackcap)
get_predCov_var_fix(fitobject,newdata=blackcap[14,],fix_X_ZAC.object=fix_X_ZAC.object)

## point predictions and variances with new X and Z
if(require("rsae", quietly = TRUE)) {
  data("landsat")
  fitobject <- HLfit(HACorn ~ PixelsCorn + PixelsSoybeans + (1|CountyName),
                     data=landsat[-33,],HLmethod="ML")
  newXandZ <- unique(data.frame(PixelsCorn=landsat$MeanPixelsCorn,
                                PixelsSoybeans=landsat$MeanPixelsSoybeans,
                                CountyName=landsat$CountyName))
  predict(fitobject,newdata=newXandZ,variances = list(predVar=TRUE))
  get_predVar(fitobject,newdata=newXandZ,variances = list(predVar=TRUE))
}

}
}
%- \keyword{ print }
