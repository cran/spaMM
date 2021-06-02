\name{predVar}
\alias{predVar}

\title{Prediction and response variances}

\description{spaMM allows computation of four variance components of prediction, returned by \code{predict} as \dQuote{...\code{Var}} attributes: \code{predVar}, \code{fixefVar}, \code{residVar}, or \code{respVar}. The phrase \dQuote{prediction variance} is used inconsistently in the literature. Often it is used to denote the uncertainty in the response (therefore, including the residual variance), but spaMM departs from this usage. Here, this uncertainly is called the response variance (\code{respVar}), while prediction variance (\code{predVar}) is used to denote the uncertainty in the linear predictor (as in Booth & Hobert, 1998). The \code{respVar} is the \code{predVar} plus the residual variance \code{residVar}.

Which components are returned is controlled in particular by the \code{type} and \code{variances} arguments of the relevant functions. \code{variances} is a list of booleans whose possible elements either match the possible returned components: \code{predVar}, \code{fixefVar}, \code{residVar}, or \code{respVar}; or may additionally include \code{linPred}, \code{disp}, \code{cov}, \code{as_tcrossfac_list} and possibly other cryptic ones. 

The \code{predict} default value for all elements is \code{NULL}, which jointly translate to no component being computed, equivalently to setting all elements to \code{FALSE}. However, setting one component to \code{TRUE} may reverse the default effect for other components. In particular, by default, component \code{predVar} implies \code{linPred=TRUE, disp=TRUE} and component \code{respVar} additionally implies \code{residVar=TRUE}; in both cases, the \code{linPred=TRUE} default by default implies \code{fixefVar=TRUE}. Calling for one variance may imply that some of its components are not only computed but also returned as a distinct attribute.    

By default the returned components are vectors of variances (with exceptions for some \code{type} value). To obtain covariance matrices (when applicable), set \code{cov=TRUE}. \code{as_tcrossfac_list=TRUE} can be used to return a list of matrices \eqn{X_i} such that the \code{predVar}\ covariance matrix equals \eqn{\sum_i X_i X'_i}. It thus provides a representation of the \code{predVar} that may be useful in particular when the \code{predVar} has large dimension, as the component \eqn{X_i}s may require less memory (being possibly non-square or sparse).  

\code{residVar=TRUE} evaluates \code{residVar} the residual variances for Gaussian or Gamma responses. 

\code{fixefVar=TRUE} evaluates \code{fixefVar}, the variance due to uncertainty in fixed effects (\bold{X}\eqn{\beta}).

Computations implying \code{linPred=TRUE} will take into account the variances of the linear predictor \eqn{\eta}, i.e. the uncertainty in fixed effects (\bold{X}\eqn{\beta}) and random effects (\bold{ZLv}), \bold{for given dispersion parameters} (see Details).
For fixed-effect models, the \code{fixefVar} calculations reduces to the \code{linPred} one.

Computations implying \code{disp=TRUE} additionally include the effect of uncertainty in estimates of dispersion parameters (\eqn{\lambda} and \eqn{\phi}), with some limitations: this effect can be computed for a scalar residual variance (\eqn{\phi}) and for several random effects with scalar variances (\eqn{\lambda}). Thus, the argument \code{variances=list(predVar=TRUE)} implies that uncertainty if linear predictor, including uncertainty in dispersion parameters, is taken into account, and the argument \code{variances=list(respVar=TRUE)} additionally includes residual variance.

}

\details{

  \code{fixefVar} is the (co)variance of \bold{X}\eqn{\beta}, deduced from the asymptotic covariance matrix of \eqn{\beta} estimates. 

  \code{linPred} is the prediction (co)variance of \eqn{\eta}=\bold{X}\eqn{\beta}+\bold{Z}\bold{v} (see \code{\link{HLfit}} Details for notation, and keep in mind that new matrices may replace the ones from the fit object when \code{newdata} are used), by default computed for given dispersion parameters. It takes into account the joint uncertainty in estimation of \eqn{\beta} and prediction of \bold{v}. 
In particular, for new levels of the random effects, \code{predVar} computation takes into account uncertainty in prediction of \bold{v} for these new levels. For \bold{prediction covariance} with a new \bold{Z}, it matters whether a single or multiple new levels are used: see Examples.

For computations implying \code{disp=TRUE}, prediction variance may also include a term accounting for uncertainty in \eqn{\phi} and \eqn{\lambda}, computed following Booth and Hobert (1998, eq. 19). This computation ignores uncertainties in spatial correlation parameters. 

\code{respVar} is the sum of \code{predVar} (pre- and post-multiplied by \eqn{\partial\mu/\partial\eta}{d\mu/d\eta} for models with non-identity link) and of \code{residVar}. 
  
These variance calculations are approximate except for LMMs, and cannot be guaranteed to give accurate results.  

}
\references{
  Booth, J.G., Hobert, J.P. (1998) Standard errors of prediction in generalized linear mixed models. J. Am. Stat. Assoc. 93: 262-272. 
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
