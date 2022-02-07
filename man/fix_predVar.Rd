\name{fix_predVar}
\alias{fix_predVar}
\title{
Prediction from models with nearly-singular covariance matrices
}
\description{
This explains how to handle a warning occurring in computation of prediction variance, where the user is directed here. 

For \bold{Matern or Cauchy} correlation models with vanishing scale factor for distances, a warning may be produced when \code{predict.HLfit} (or \code{get_predVar}, etc.) is called with non-NULL \code{newdata}, because a nearly-singular correlation matrix of the random effect is met. \bold{To decide what to do} in that case, users should compare the values of \code{get_predVar(.)} and \code{get_predVar(., newdata=myfit$data)}(see Example below). In the absence of numerical inaccuracies, The two values should be identical, and in the presence of such inaccuracies, the more reliable value is the first one. In really poor cases, the second syntax may yield negative prediction variances. If users deem the inaccuracies too large, they should use \code{control=list(fix_predVar=TRUE)} in the next call to \code{predict.HLfit} (or \code{get_predVar}, etc.) as shown in the Example. The drawback of this control is that the computation may be slower, and might even exceed memory capacity for large problems (some matrix operations being performed with exact rational arithmetic, which is memory-consuming for large matrices). it is also still experimental, in the sense that I fear that bugs (\code{stop}) may occur. If the user instead chooses \code{control=list(fix_predVar=FALSE)}, the default standard floating-point arithmetic is used, but no warning is issued. 

For \code{fix_predVar} left NULL (the default), standard floating-point arithmetic is also used. But in addition (with exceptions: see Details), the warning keeps being issued, and the (possibly costly) computation of the inverse of the correlation matrix is not stored in the fitted model object, hence is repeated for each new prediction variance computation. This is useful to remind users that something needs to be done, but for programming purposes where repeated warnings may be a nuisance, one can use \code{control=list(fix_predVar=NA)} which will issue a warning then perform as \code{control=list(fix_predVar=FALSE)}, i.e. store an approximate inverse so the warning is not issued again. Finally, \code{control=list(fix_predVar=NaN)} will remove the inverse of the correlation matrix from the fitted model object, and start afresh as if the control was NULL.

}

\details{

Nearly-singular correlation matrices of random effects occur in several contexts. For random-slope models, it commonly occurs that the fitted correlation between the random effects for Intercept and slope is 1 or -1, in which case the correlation matrix between these random effects is singular. This led to quite inaccurate computations of prediction variances in spaMM prior to version 3.1.0, but this problem has been fixed. 

\code{control=list(fix_predVar=NaN)} may be more appropriate than \code{control=list(fix_predVar=NULL)} when \code{predict.HLfit} is called through code that one cannot control. For this reason, spaMM provides another mode of control of the default. It will convert \code{control=list(fix_predVar=NULL)} to other values when the call stack has call names matching the patterns given by \code{spaMM.getOption("fix_predVar")} (as understood by \code{\link{grep}}). Thus if \code{spaMM.getOption("fix_predVar")$"NA"=="MSL|bboptim"}, the default behaviour is that defined by \code{control=list(fix_predVar=NA)} when \code{predict.HLfit} is called through \code{Infusion::MSL} or \code{blackbox::bboptim}. FALSE or TRUE are handled in a similar way.

}
\examples{
data("blackcap")
fitobject <- corrHLfit(migStatus ~ 1 + Matern(1|longitude+latitude),data=blackcap,
                       ranFix=list(nu=10,rho=0.001)) ## numerically singular C
get_predVar(fitobject,newdata=blackcap[6,]) 
## => warning => let us apply the recommended procedure:
get_predVar(fitobject) 
get_predVar(fitobject,newdata=fitobject$data) 
# Negative values again in the second case => easy decision:
get_predVar(fitobject,newdata=blackcap[1:6,], 
            control=list(fix_predVar=TRUE)) # now it's accurate
            # and the accuracy control is stored in the object:
get_predVar(fitobject,newdata=blackcap[1:6,]) 
# Clean and start afresh:
get_predVar(fitobject,newdata=blackcap[1:6,], 
            control=list(fix_predVar=NaN)) 
}
