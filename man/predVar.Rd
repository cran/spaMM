\name{predVar}
\alias{predVar}

\title{Prediction and response variances}

\description{spaMM allows computation of four variance components of prediction, returned by \code{predict} as \dQuote{...\code{Var}} attributes: \code{predVar}, \code{fixefVar}, \code{residVar}, or \code{respVar}. The phrase \dQuote{prediction variance} is used inconsistently in the literature. Often it is used to denote the uncertainty in the response (therefore, including the residual variance), but \pkg{spaMM} follows some literature for mixed models in departing from this usage. Here, this uncertainly is called the response variance (\code{respVar}), while prediction variance (\code{predVar}) is used to denote the uncertainty in the linear predictor (as in Booth & Hobert, 1998; see also Jeske & Harville, 1988). The \code{respVar} is the \code{predVar} plus the residual variance \code{residVar}.

Which components are returned is controlled in particular by the \code{type} and \code{variances} arguments of the relevant functions. \code{variances} is a list of booleans whose possible elements either match the possible returned components: \code{predVar}, \code{fixefVar}, \code{residVar}, or \code{respVar}; or may additionally include \code{linPred}, \code{disp}, \code{cov}, \code{as_tcrossfac_list} and possibly other cryptic ones. 

The \code{predict} default value for all elements is \code{NULL}, which jointly translate to no component being computed, equivalently to setting all elements to \code{FALSE}. However, setting one component to \code{TRUE} may reverse the default effect for other components. In particular, by default, component \code{predVar} implies \code{linPred=TRUE, disp=TRUE} and component \code{respVar} additionally implies \code{residVar=TRUE}; in both cases, the \code{linPred=TRUE} default by default implies \code{fixefVar=TRUE}. Calling for one variance may imply that some of its components are not only computed but also returned as a distinct attribute.    

By default the returned components are vectors of variances (with exceptions for some \code{type} value). To obtain covariance matrices (when applicable), set \code{cov=TRUE}. \code{as_tcrossfac_list=TRUE} can be used to return a list of matrices \eqn{X_i} such that the \code{predVar} covariance matrix equals \eqn{\sum_i X_i X'_i}. It thus provides a representation of the \code{predVar} that may be useful in particular when the \code{predVar} has large dimension, as the component \eqn{X_i}s may require less memory (being possibly non-square or sparse).  

\code{residVar=TRUE} evaluates \code{residVar} the residual variance. For families without a dispersion parameter (e.g., \code{binomial} or \code{poisson}), this is as given by the \code{variance} function of the \code{family} object (in the binomial case, it is thus the variance of a single binary draw). For families with a dispersion parameter (such as \eqn{\phi} for gaussian or Gamma families, negative-binomial, beta), it is the residual variance as function of the dispersion parameter, whether this parameter is a single scalar or follows a more complex residual-dispersion model.  Prior weights are however ignored (see the \code{\link{residVar}} etractor for the opposite feature). For the beta-binomial family, it is also the variance of a single binary draw; although this family has a residual-dispersion parameter the latter variance is not affected by it.

\code{fixefVar=TRUE} evaluates \code{fixefVar}, the variance due to uncertainty in fixed effects (\bold{X}\eqn{\beta}).

Computations implying \code{linPred=TRUE} will take into account the variances of the linear predictor \eqn{\eta}, i.e. the uncertainty in fixed effects (\bold{X}\eqn{\beta}) and random effects (\bold{ZLv}), \bold{for given dispersion parameters} (see Details).
For fixed-effect models, the \code{fixefVar} calculations reduces to the \code{linPred} one.

Computations implying \code{disp=TRUE} additionally include the effect of uncertainty in estimates of dispersion parameters (\eqn{\lambda} and \eqn{\phi}), with some limitations (see Details). \code{variances=list(predVar=TRUE)}, which evaluates the uncertainty of linear predictor, implies \code{disp=TRUE} by default, meaning that it includes such effects of uncertainty in dispersion parameters on the linear predictor. \code{variances=list(respVar=TRUE)} performs similarly but additionally includes the residual variance in the returned variance.
}

\details{

  \code{fixefVar} is the (co)variance of \bold{X}\eqn{\beta}, deduced from the asymptotic covariance matrix of \eqn{\beta} estimates. 

  \code{linPred} is the prediction (co)variance of \eqn{\eta}=\bold{X}\eqn{\beta}+\bold{Z}\bold{v} (see \code{\link{HLfit}} Details for notation, and keep in mind that new matrices may replace the ones from the fit object when \code{newdata} are used), by default computed for given dispersion parameters. It takes into account the joint uncertainty in estimation of \eqn{\beta} and prediction of \bold{v}. 
In particular, for new levels of the random effects, \code{predVar} computation takes into account uncertainty in prediction of \bold{v} for these new levels. For \bold{prediction covariance} with a new \bold{Z}, it matters whether a single or multiple new levels are used: see Examples.

For computations implying \code{disp=TRUE}, prediction variance may also include a term accounting for uncertainty in \eqn{\phi} and \eqn{\lambda}, computed following Booth and Hobert (1998, eq. 19). This computation acheives its originally described purpose for a scalar residual variance (\eqn{\phi}) and for several random effects with scalar variances (\eqn{\lambda}). This computation ignores uncertainties in spatial correlation parameters. 

The (1998) formulas for the effect of uncertainty in dispersion parameters are here also applied to the variance parameters in random-coefficient terms, but with a one-time warning. Not only this is not expected to account for the uncertainty of the correlation parameter(s) of such terms, but the result will depend on the internal representation (the \dQuote{square root}) of the covariance matrix, which may differ among the different fitting algorithms that may be used by \pkg{spaMM}. % ___F I X M E____ TAG rc_dispcov 

\code{respVar} is the sum of \code{predVar} (pre- and post-multiplied by \eqn{\partial\mu/\partial\eta}{d\mu/d\eta} for models with non-identity link) and of \code{residVar}. 
  
These variance calculations are approximate except for LMMs, and cannot be guaranteed to give accurate results.  

}
\references{
  Booth, J.G., Hobert, J.P. (1998) Standard errors of prediction in generalized linear mixed models. J. Am. Stat. Assoc. 93: 262-272. 
  
  Jeske, Daniel R. & Harville, David A. (1988) Prediction-interval procedures and (fixed-effects) confidence-interval procedures for mixed linear models. Communications in Statistics - Theory and Methods, 17: 1053-1087. \doi{10.1080/03610928808829672} 
}
\examples{
\dontrun{
# (but run in help("get_predVar"))
data("blackcap")
fitobject <- fitme(migStatus ~ 1 + Matern(1|longitude+latitude),data=blackcap,
                       fixed=list(nu=4,rho=0.4,phi=0.05))

#### multiple controls of prediction variances
# (1) fit with an additional random effect
grouped <- cbind(blackcap,grp=c(rep(1,7),rep(2,7))) 
fitobject <- fitme(migStatus ~ 1 +  (1|grp) +Matern(1|longitude+latitude),
                       data=grouped,  fixed=list(nu=4,rho=0.4,phi=0.05))

# (2) re.form usage to remove a random effect from point prediction and variances: 
predict(fitobject,re.form= ~ 1 +  Matern(1|longitude+latitude))

# (3) comparison of covariance matrices for two types of new data
moregroups <- grouped[1:5,]
rownames(moregroups) <- paste0("newloc",1:5)
moregroups$grp <- rep(3,5) ## all new data belong to an unobserved third group 
cov1 <- get_predVar(fitobject,newdata=moregroups,
                     variances=list(linPred=TRUE,cov=TRUE))
moregroups$grp <- 3:7 ## all new data belong to distinct unobserved groups
cov2 <- get_predVar(fitobject,newdata=moregroups,
                     variances=list(linPred=TRUE,cov=TRUE))
cov1-cov2 ## the expected off-diagonal covariance due to the common group in the first fit.

}
## see help("get_predVar") for further examples
}
