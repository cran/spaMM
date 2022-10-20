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
  Prediction from a model fit
}
\description{
The following functions can be used to compute point predictions and/or various measures of uncertainty associated to such predictions. \code{predict} can be used for prediction of the response variable by its expected value obtained as (the inverse link transformation of) the linear predictor (\eqn{\eta}) and more generally for terms of the form \bold{X}_n\eqn{\beta}+\bold{Z}_n\bold{L}\bold{v}, for new design matrices \bold{X}_n and \bold{Z}_n. Various components of prediction variances and predictions intervals can also be computed using \code{predict}. The \code{get_}... functions are convenient extractors for such components. \code{get_predCov_var_fix} extracts a block of a prediction covariance matrix. It was conceived for the specific purpose of computing the spatial prediction covariances between two \dQuote{new} sets of geographic locations, without computing the full covariance matrix for both the new locations and the original (fitted) locations. When one of the two sets of new locations is fixed while the other varies, some expensive computations can be performed once for all sets of new locations, and be provided as the \code{fix_X_ZAC.object} argument. The \code{preprocess_fix_corr} extractor is designed to compute this argument. 
}
\usage{
\method{predict}{HLfit}(object, newdata = newX, newX = NULL, re.form = NULL,
                        variances=list(), binding = FALSE, intervals = NULL,
                        level = 0.95, blockSize = 2000L, type = "response", 
                        verbose=c(showpbar=eval(spaMM.getOption("barstyle"))), 
                        control=list(), ...)
get_predCov_var_fix(object, newdata = NULL, fix_X_ZAC.object, fixdata, re.form = NULL,
                    variances=list(disp=TRUE,residVar=FALSE,cov=FALSE), 
                    control=list(),  ...)    
preprocess_fix_corr(object, fixdata, re.form = NULL,
                   variances=list(residVar=FALSE, cov=FALSE), control=list())
get_fixefVar(...)
get_predVar(..., variances=list(), which="predVar")
get_residVar(...)
get_respVar(...)
get_intervals(..., intervals="predVar")
}
\arguments{
  \item{object}{
  The return object of fitting functions \code{HLfit,corrHLfit,HLCor}... returning an object inheriting from \code{HLfit} class.
}
  \item{newdata}{
  \bold{Either} NULL, a matrix or data frame, or a numeric vector. 
  
  If \code{NULL}, the original data are reused. Otherwise, all variables required to evaluate model formulas must be included. Which variables are required may depend on other arguments: see \dQuote{prediction with given phi's} example, also illustrating the syntax when formulas include an offset.  
  
  If \code{newdata} is a numeric vector, its names (if any) are ignored. This makes it easier to use \code{predict} as an objective function for an 
  optimization procedure such as \code{optim}, which calls the objective function on unnamed vectors. However, one must make sure that the order of elements in the vector is the order of first occurrence of the variables in the model formula. This order can be checked in the error message returned when calling \code{predict} on a \code{newX} vector of clearly wrong size, e.g. \code{predict(<object>,newdata=numeric(0))}.
}
 \item{newX}{equivalent to newdata, available for back-compatibility} 
 \item{re.form}{
  formula for random effects to include.  By default, it is NULL, in which case all random effects are included. If it is NA, no random effect is included. If it is a formula, only the random effects it contains are retained. The other variance components are removed from both point prediction and \code{variances} calculations. If you want to retain only the spatial effects in the point prediction, but all variances, either use re.form  and add missing variances (on linear predictor scale) manually, or ignore this argument and see Details and Examples for different ways of controlling variances.     
 }
 \item{variances}{
   A list whose elements control the computation of different estimated variances. 
   \code{predict} can return four components of prediction variance: \code{fixefVar}, \code{predVar}, \code{residVar} and \code{respVar}. 
   They are all returned as attributes of the point predictions. 
   
   In particular, \code{variances=list(predVar=TRUE)} is suitable for uncertainty in point prediction, distinguished from the response variance given by \code{list(respVar=TRUE)}. See the \code{\link{predVar}} help page for further explanations and other options. 
 }
 \item{intervals}{
   NULL or character string or vector of strings. Provides prediction intervals with nominal level \code{level}, deduced from the given prediction variance term, e.g. \code{intervals="predVar"}. Currently only intervals from \code{fixefVar} and \code{predVar} (and for LMMs \code{respVar} including the residual variance) may have a probabilistic meaning. Intervals returned in other cases are (currently) meaningless. 
 }
 \item{which}{any of \code{"predVar"},\code{"respVar"},\code{"residVar"}, \code{"fixefVar"}, \code{"intervals"}, or \code{"naive"} }
 \item{level}{Coverage of the intervals.}
 \item{binding}{
  If \code{binding} is a character string, the predicted values are bound with the \code{newdata} and the result is returned as a data frame. The predicted values column name is the given \code{binding}, or a name based on it if the \code{newdata} already include a variable with this name. 
  If \code{binding} is \code{FALSE}, The predicted values are returned as a one-column matrix and the data frame used for prediction is returned as an attribute
  (unless it was \code{NULL}). If \code{binding} is \code{NA}, a vector is returned, without the previous attributes.
 }
 \item{fixdata}{A data frame describing reference data whose covariances with variable \code{newdata} may be requested.}
 \item{fix_X_ZAC.object}{The return value of calling \code{preprocess_fix_corr} (see trivial Example). This is a more efficient way of providing information about the \code{fixdata} for repeated calls to \code{get_predCov_var_fix} with variable \code{newdata}.}
 \item{blockSize}{
   Mainly for development purposes. For original or new data with many rows, it may be more efficient to split these data in small blocks, and this gives the maximum number or rows of the blocks. However, this will be ignored if a prediction covariance matrix is requested.    
 }
\item{type}{character string; The returned point predictions are on the response scale if \code{type="response"} (the default; for binomial response, a frequency 0<.<1). It is on the linear predictor scale if \code{type="link"}. \cr
* The \dQuote{prediction variance} (as opposed to the response variance, see \code{\link{predVar}}) that may be returned as a \code{"predVar"}  attribute of the point predictions is always on the linear predictor scale, even when \code{type="response"}. If you want to extract this \code{predVar} transformed to the response scale, use \code{predict(.,variances=list(respVar=TRUE))} and take the difference between the \code{respVar} and \code{residVar} attributes of the result.\cr
* Prediction intervals (as opposed to the response intervals) will be on the linear predictor or response scale depending on \code{type} (new to versions more recent than 3.12.0).}
\item{control}{A list; a warning will direct you to relevant usage when needed.}
\item{verbose}{A vector of booleans; it single currently used element is \code{"showpbar"}, which controls whether to show a progress bar in certain prediction variance computations.}
\item{\dots}{
  further arguments passed to or from other methods. For the \code{get_}... functions, they are passed to \code{predict}. 
}
}
\value{
  See Details in \code{\link{Tpoisson}} for questions specific to truncated distributions.
  
  For \code{predict}, a matrix or data frame (according to the \code{binding} argument), with optional attributes \code{frame}, \code{intervals}, \code{predVar}, \code{fixefVar}, \code{residVar}, and/or \code{respVar}, the last four holding one or more variance vector or covariance matrices. The further attribute \code{fittedName} contains the binding name, if any.  
  
  The \code{get_}... extractor functions call \code{predict} and extract from its result the attribute implied by the name of the extractor. By default, \code{get_intervals} will return prediction intervals using \code{predVar}.
  \code{get_predVar} with non-default \code{which} argument has the same effect as the \code{get_}... function whose name is implied by \code{which}.
}
\details{
  See the \code{\link{predVar}} help page for information about the different concepts of prediction variances handled by spaMM (uncertainty of point prediction vs. of response) and about options controlling their computation.

  If  \code{newdata} is NULL, \code{predict} returns the fitted responses, including random effects, from the object. 
  Otherwise it computes new predictions including random effects as far as possible.   
  For spatial random effects it constructs a correlation matrix \bold{C} between new locations and locations in the original fit. Then it infers the random effects in the new locations as   \bold{C} (\bold{L}'\eqn{)^{-1}} \bold{v} (see \code{\link{spaMM}} for notation).   For non-spatial random effects, it checks whether any group (i.e., level of a random effect) in the new data was represented in the original data, and it adds the inferred random effect for this group to the prediction for individuals in this group. 

In the \bold{point prediction} of the linear predictor, the unconditional expected value of \eqn{u} is assigned to the realizations of \eqn{u} for unobserved levels of non-spatial random effects (it is zero in GLMMs but not for non-gaussian random effects), and the inferred value of \eqn{u} is assigned in all other cases. Corresponding values of \eqn{v} are then deduced. This computation yields the classical \dQuote{BLUP} or empirical Bayes predictor in LMMs, but otherwise it may yield less well characterized predictors, where \dQuote{unconditional} \eqn{v} may not be its expected value when the \code{rand.family} link is not identity. 

\bold{Intervals} computations use the relevant variance estimates plugged in a Gaussian approximation, except for the simple linear model where it uses Student's \emph{t} distribution.   

}
\seealso{
\code{\link{predVar}} for information specific to prediction variances sensu lato;
\code{\link{get_cPredVar}} for a bootstrap-corrected version of \code{get_predVar};
\code{\link{residVar}} for an alternative extractor for residual variances, more general than \code{get_residVar}.
}

\examples{
data("blackcap")
fitobject <- fitme(migStatus ~ 1 + Matern(1|longitude+latitude),data=blackcap,
                       fixed=list(nu=4,rho=0.4,phi=0.05))
predict(fitobject)

#### multiple controls of prediction variances
## (1) fit with an additional random effect
grouped <- cbind(blackcap,grp=c(rep(1,7),rep(2,7))) 
fitobject2 <- fitme(migStatus ~ 1 +  (1|grp) +Matern(1|longitude+latitude),
                       data=grouped,  fixed=list(nu=4,rho=0.4,phi=0.05))

## (2) re.form usage to remove a random effect from point prediction and variances: 
predict(fitobject2,re.form= ~ 1 +  Matern(1|longitude+latitude))

## (3) comparison of covariance matrices for two types of new data
moregroups <- grouped[1:5,]
rownames(moregroups) <- paste0("newloc",1:5)
moregroups$grp <- rep(3,5) ## all new data belong to an unobserved third group 
cov1 <- get_predVar(fitobject2,newdata=moregroups,
                     variances=list(linPred=TRUE,cov=TRUE))
moregroups$grp <- 3:7 ## all new data belong to distinct unobserved groups
cov2 <- get_predVar(fitobject2,newdata=moregroups,
                     variances=list(linPred=TRUE,cov=TRUE))
cov1-cov2 ## the expected off-diagonal covariance due to the common group in the first fit.

\dontrun{
#### Other extractors:
#
fix_X_ZAC.object <- preprocess_fix_corr(fitobject,fixdata=blackcap)
#
# ... for use in multiple calls to get_predCov_var_fix():
#
get_predCov_var_fix(fitobject,newdata=blackcap[14,],fix_X_ZAC.object=fix_X_ZAC.object)

#### Prediction with distinct given phi's in different locations, 
#   as specified by a resid.model:
#
varphi <- cbind(blackcap,logphi=runif(14))
vphifit <- fitme(migStatus ~ 1 + Matern(1|longitude+latitude), 
                     resid.model = list(formula=~0+offset(logphi)),
                     data=varphi,  fixed=list(nu=4,rho=0.4))
#
# For respVar computation (i.e., response variance, often called prediction variance), 
#   one then also needs to provide the variables used in 'resid.model', here 'logphi':
#
get_respVar(vphifit,newdata=data.frame(latitude=1,longitude=1,logphi=1))
#
# For default 'predVar' computation (i.e., uncertainty in point prediction), 
#   this is not needed:
#
get_predVar(vphifit,newdata=data.frame(latitude=1,longitude=1))                     

#### point predictions and variances with new X and Z
#
if(requireNamespace("rsae", quietly = TRUE)) {
  data("landsat", package = "rsae")
  fitobject <- fitme(HACorn ~ PixelsCorn + PixelsSoybeans + (1|CountyName),
                     data=landsat[-33,])
  newXandZ <- unique(data.frame(PixelsCorn=landsat$MeanPixelsCorn,
                                PixelsSoybeans=landsat$MeanPixelsSoybeans,
                                CountyName=landsat$CountyName))
  predict(fitobject,newdata=newXandZ,variances = list(predVar=TRUE))
  get_predVar(fitobject,newdata=newXandZ,variances = list(predVar=TRUE))
}

}
}
%- \keyword{ print }
