\name{corrFamily}
\alias{corrFamily}
\title{
Using corrFamily constructors and descriptors.
}
\description{
One can declare and fit correlated random effects belonging to a user-defined correlation (or covariance) model (i.e., a parametric family of correlation matrices, although degenerate case with no parameter are also possible). This documentation is a first introduction to this feature. It is experimental in the sense that its design has been formalized only from a limited number of corrFamily examples, and that the documentation is not mature. Implementing prediction for random-effects defined in this way may be tricky. A distinct documentation \code{\link{corrFamily-design}} provides more information for the efficient design of new correlation models to be fitted in this way.  

A simple example of random-effect model implemented in this way is the autoregressive model of order \eqn{p} (AR(p) in the literature; specifically documented elsewhere, see \code{\link{ARp}}). It can be used as a formula term like other autocorrelated random-effects predefined in \pkg{spaMM}, to be fitted by \code{fitme}
or \code{fitmv}:
\preformatted{
fitme(lh ~ 1 + ARp(1|time, p=3),  # <= declaration of random effect
  < data and other possible arguments >)
}
User-defined correlation models should be registered for this simple syntax to work (see Details for an alternative syntax):  
\preformatted{
myARp <- ARp                   # 'myARP' is thus a user-defined model
register_cF("myARp")        # Register it so that the next call works
fitme(lh ~ 1 + myARp(1|time, p=3),  
  < data and other possible arguments >)
}
The \code{ARp} object here copied in \code{myARp} is a function (the \emph{corrFamily constructor}) which returns a \code{list} (the \emph{corrFamily descriptor}) which contains the necessary information to fit a random effect with an AR(p) correlation. The \code{p} argument in the \code{myARp(1|time, p=3)} term enforces evaluation of \code{myARp(p=3)}, producing the descriptor for the AR(3) model. The structure of this descriptor is 
\preformatted{
List of 5
 $ Cf            :function (parvec)  
  ..- < with some attribute >
 $ tpar         : num [1:3] 0.5 0.333 0.25
 $ type         : chr "precision"
 $ initialize   :function (Zmatrix, ...) 
  ..- < with some attribute >
 $ fixed        : NULL
 $ calc_moreargs:function (corrfamily, ...)  
  ..- < with some attribute >
 $ levels_type        : chr "time_series"
 $ calc_corr_from_dist:function (ranFix, char_rd, distmat, ...)  
  ..- < with some attribute >
 < and possibly other elements > 
}
The meaning of these elements and some additional ones is explained below. 

Only \code{Cf} and \code{tpar} are necessary elements of a corrFamily object. If one designs a new descriptor where some other elements are absent, \pkg{spaMM} will try to provide plausible defaults for these elements. Further, if the descriptor does not provide parameter names (as the names of \code{tpar}, or in some more cryptic way), default names \code{"p1"}, \code{"p2"}... will be provided.    
}

\usage{
## corrFamily descriptor provided as a list of the form
#
# list(Cf=<.>, tpar=<.>, ...)

## corrFamily constructor: any function that returns 
#    a valid corrFamily descriptor
#
# function(tpar=<.>, fixed=<.>, ...) # typical but not mandatory arguments

## There is a distinct documentation page for 'register_cF'.
}
\arguments{

\bold{\emph{Elements} of the corrFamily \emph{descriptor}}:

\item{Cf}{(required): function returning the correlation matrix (or covariance matrix, or their inverse), given its first argument, a parameter vector.}
\item{tpar}{(required): a feasible argument of \code{Cf}. \code{tpar} is \bold{not} an initial \bold{nor} a fixed value.}
\item{type}{optional, but required if the return value of \code{Cf} is an inverse correlation matrix rather than a correlation matrix, in which case one should specify \code{type="precision"}.}
\item{fixed}{optional: fixed values for some correlation parameters, provided as a named vector with names consistent with those to be used for \code{tpar}. This is conceived to achieve the same statistical fit as by using the \code{fixed} argument of \code{fitme}, although the structure of the result of the fit differs in some subtle ways whether parameters are fixed through the descriptor or through the fitting function (see Examples in \code{\link{ARp}}).}
\item{calc_moreargs}{optional: a function returning a list with possible elements \code{init}, \code{lower} and \code{upper} for control of estimation (and possibly other elements for other purposes). If the descriptor does not provide this function, a default \code{calc_moreargs} will be provided, implementing unbounded optimization.}
\item{initialize}{optional: a function evaluating variables that may be needed repeatedly by \code{Cf} or \code{Af}.}
\item{Af}{This function should be defined if the correlation model requires an \bold{A} matrix (the middle term in the case the design matrix of a random effect term is described by a triple matrix product \bold{ZAL} as described in \code{\link{random-effects}}). Examples can be found in the descriptors returned by the \code{\link{ranGCA}} and \code{\link{MaternIMRFa}} constructors. 
}
\item{levels_type}{In the above example its value \code{"time_series"} informs \pkg{spaMM} that levels of the random effect should be considered for all integer values within the range of the \code{time} variable, not only for levels present in the data. If this element is not provided by the constructor, \pkg{spaMM} will internally assume a \code{levels_type} suitable for geostatistical models. Further level types may be defined in the future.}
\item{calc_corr_from_dist, make_new_corr_lists}{Functions possibly needed for prediction (see Details).} 
\item{need_Cnn}{optional: a boolean; default is TRUE. Controls prediction computations (see Details).}
\item{public}{An environment where some variables can be saved, typically by the \code{initialize} expression, for inspection at user level and for re-use. See \code{\link{diallel}} for an example.}

\bold{Optional \emph{fitting-function} arguments}:

\item{lower,upper,init}{(possibly incomplete) parameter vectors to control parameter estimation. Parameter names (consistent with those to be used for the \code{tpar} argument) may be required to disambiguate incomplete vectors. Any parameter value not specified in this way will be sought in the return value of the \code{calc_moreargs} function, if provided in the descriptor.} 
%(given the distinctive use of \code{tpar}, \code{spaMM} never tries to use it as initial value).
\item{fixed}{(optional): vector of fixed values. Parameter names may be needed as for \code{lower,upper,init}.}

\bold{\emph{Arguments} of the corrFamily \emph{constructor}}

These may be ad libitum, as design rules are defined only for the returned descriptor. However, arguments \code{tpar}, \code{fixed}, and \code{public}  predefined constructors, such as \code{ARp}, are designed to match the above respective elements of the descriptor. 
}
\details{
\bold{* Constructor elements for prediction:}

For prediction of autocorrelated random effects, one must first assess whether levels of the random effect not represented in the fit are possible in new data (corresponding to new spatial locations in geostatistical models, or new time steps in time series). In that case \code{need_Cnn} must be TRUE (Interpolated MRFs do not require this as all required random-effect levels are determined by the IMRF \code{mesh} argument rather than by the fitted data or new data).

Next, for autocorrelated random effects where \code{need_Cnn} is TRUE, a \code{make_new_corr_lists} function must be provided, except when a \code{calc_corr_from_dist} function is instead provided (which may be sufficient for models that construct the correlation from a spatial distance matrix). When \code{need_Cnn} is FALSE, a \code{make_new_corr_lists} function may still be needed. 

The Examples section provides a simple example of such design, and the source code of the \code{ARp} or \code{ARMA} constructors provide further examples. They show that the \code{make_new_corr_lists} function may assign matrices or vectors as elements of three lists contained in a \code{newLv_env} environment. A matrix is assigned in the \code{cov_newLv_oldv_list}, specifying correlations between \dQuote{new} levels of the random effect (implied by the new data) and \dQuote{old} levels (those already included in the design matrix of the random effect for the fit). If \code{need_Cnn} is TRUE, a second matrix may be assigned in the \code{cov_newLv_newLv_list}, specifying correlation between \dQuote{new} levels, and the diagonal of this matrix is assigned in the \code{diag_cov_newLv_newLv_list}. The overall structure of the code (the conditions where these assignments are performed, and the list indices), should be conserved. 


\bold{* Fitting a \code{corrFamily} without a constructor:}

It is possible to use an unregistered corrFamily, as follows: 
\preformatted{
AR3 <- ARp(p=3)          # Generate descriptor of AR model of order 3

fitme(lh ~ 1 + corrFamily(1|time),  # <= declaration of random effect
  covStruct=list( 
    corrFamily= AR3     # <= declaration of its correlation structure
  ), 
  < data and other possible arguments >)
}

Here the fit only uses a descriptor list, not a constructor function. This descriptor is here provided to the fitting function as an element of the \code{\link{covStruct}} argument (using the general syntax of this argument), and in the model formula the corresponding random effect is declared as a term of the form\cr
\code{corrFamily(1|<grouping factor>)}. 

This syntax is more complex than the one using a registered constructor, but it might be useful for development purposes (one only has to code the descriptor, not the constructor function). However, it is not general; in particular, using registered constructors may be required to obtain correct results when fitting multivariate-response models by \code{fitmv}.
}
\seealso{
See \code{\link{ARp}}, \code{\link{diallel}}, and \code{\link{MaternIMRFa}} for basic examples of using a predefined corrFamily descriptor, and \code{\link{corrFamily-design}} for more geeky stuff including examples of implementing simple new correlation families. 
}
\examples{
\dontrun{
### Minimal (with many features missing) reimplementation 
#     of corrMatrix() terms as a corrFamily 


corrMatrix_cF <- function(corrMatrix) {
  
  force(corrMatrix) # Makes it available in the environment of the functions next defined.
  oldZlevels <- NULL
  
  initialize <- function(Zmatrix, ...) {
    oldZlevels <<- colnames(Zmatrix) # Pass info about levels of the random effect in the data.
  }
  
  Cf <- function(newlevels=oldZlevels ) {
    if (length(newlevels)) {
      corrMatrix[newlevels,newlevels]
    } else corrMatrix[oldZlevels,oldZlevels] # for Cf(tpar=numeric(0L))
  }
  
  calc_moreargs <- function(corrfamily, ...) {
    list(init=c(),lower=c(),upper=c())
  }
  
  make_new_corr_lists <- function(newLv_env, which_mats, newZAlist, new_rd, ...) {
    newlevels <- colnames(newZAlist[[new_rd]])
    newLv_env$cov_newLv_oldv_list[[new_rd]] <- corrMatrix[newlevels,oldZlevels, drop=FALSE]
    if (which_mats$nn[new_rd]) {
      newLv_env$cov_newLv_newLv_list[[new_rd]] <- corrMatrix[newlevels,newlevels, drop=FALSE]
    } else { 
      newLv_env$diag_cov_newLv_newLv_list[[new_rd]] <- rep(1,length(newlevels)) 
    }
  }
  
  list(Cf=Cf, tpar=numeric(0L), initialize=initialize, calc_moreargs=calc_moreargs, 
       make_new_corr_lists=make_new_corr_lists,
       tag="corrMatrix_cF") 
}

register_cF("corrMatrix_cF")

# usage:

data("blackcap") 
MLcorMat <- MaternCorr(proxy::dist(blackcap[,c("latitude","longitude")]),
                       nu=0.6285603,rho=0.0544659)
corrmat <- proxy::as.matrix(MLcorMat, diag=1)

fitme(migStatus ~ means+ corrMatrix_cF(1|name, corrMatrix=corrmat),data=blackcap,
      corrMatrix=MLcorMat,method="ML")
      
unregister_cF("corrMatrix_cF") # Tidy things before leaving.         

}
}
