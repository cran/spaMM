\name{get_cPredVar}
\alias{get_cPredVar}
\title{
Estimation of prediction variance with bootstrap correction
}
\description{
This function is similar to \code{\link{get_predVar}} except that is uses a bootstrap procedure to correct for bias in the evaluation of the prediction variance.     
}
\usage{
get_cPredVar(pred_object, newdata, nsim, seed, type = "residual", 
             variances=NULL, nb_cores = NULL, fit_env = NULL,
             sim_object=pred_object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pred_object}{an object of class \code{HLfit}, as returned by the fitting functions in \code{spaMM}.}
  \item{newdata}{passed to \code{\link{predict.HLfit}} (it thus represents a prediction design, not to be confused with the bootstrap samples)}
  \item{nsim}{passed to \code{\link{simulate.HLfit}}}
  \item{seed}{passed to \code{\link{simulate.HLfit}}}
  \item{type}{passed to \code{\link{simulate.HLfit}}}
  \item{variances}{NULL or list; \code{variances["cov"]} will be passed to \code{\link{predict.HLfit}} to control whether a covariance matrix is computed or not. Other elements are currently ignored.}
  \item{nb_cores}{integer: number of cores to use for parallel computation of bootstrap. The default is \code{spaMM.getOption("nb_cores")}, and 1 if the latter is NULL. \code{nb_cores=1} prevents the use of parallelisation procedures.}
  \item{fit_env}{For parallel computations: an environment containing objects to be passed to the cores. They should have the same name in \code{fit_env} as in the environment they are passed from.}
  \item{sim_object}{an object of class \code{HLfit}, passed to \code{\link{simulate.HLfit}} as its \code{object} argument. Simulating from this object must produce response values that can be used as replacement to those of the original fitted \code{pred_object}. In standard usage, \code{sim_object=pred_object} (the default).}
}
\details{
The result provided by \code{get_cPredVar} is similar to the CMSEP (Conditional Mean Standard Error of Prediction) introduced by Booth and Hobert (1998; \dQuote{B&H}). This paper is known for pointing the importance of using conditional variances when they differ from unconditional ones. This is hard to miss in spatial models, where the relevant prediction variance typically depends on the variance of random effects conditional on the data. Thus, the alternative function \code{get_predVar} already accounts for this and returns a prediction variance that depends on a joint covariance of fixed-effect estimates and of random effects given the data. 

B&H also used a conditional bootstrap procedure to correct for some bias. \code{get_cPredVar} implements a similar procedure, in contrast to \code{get_predVar}. Their conditional bootstrap procedure is not applicable for autocorrelated random effects, and parametric bootstrapping of the residuals of the fitted model (as implied by the default value of argument \code{type}) is used instead here. Apart from this difference, the returned value includes exactly the same terms as those discussed by B&H: their \dQuote{naive estimate} \eqn{\nu_i} and its bootstrap correction \eqn{b_i}, their correction \eqn{\beta} for uncertainty in fixed-effect coefficients, and their correction \eqn{\sigma^2} for uncertainty in dispersion parameters. 

The \dQuote{naive estimate} \eqn{\nu_i} is not generally an estimate of anything uniquely defined by the model parameters: for correlated random effects, it depends on the \dQuote{root} of the correlation matrix of the random effects, which is not unique. Thus \eqn{\nu_i} is not unique, and may differ for example for equivalent fits by sparse-precision methods vs. other methods. Nevertheless, \code{attr(cpredvar,"info")$naive} does recover published values in the Examples below, as they involve no correlation matrix. 
}
\value{
A vector of prediction variances, with an attribute \code{info} which is an \bold{environment} containing variables:
\item{SEs}{the standard errors of the estimates (which are those of the bootstrap replicates)}
\item{bias}{the bias term}
\item{maive}{B&H's \dQuote{naive} \eqn{\nu_i}}
}
\references{
  Booth, J.G., Hobert, J.P. (1998) Standard errors of prediction in generalized linear mixed models. J. Am. Stat. Assoc. 93: 262-272. 
}
\examples{
\dontrun{
if(requireNamespace("rsae", quietly = TRUE)) {
  # LMM example from Booth & Hobert 1998 JASA
  data("landsat", package = "rsae")
  fitCorn <- fitme(HACorn ~ PixelsCorn + PixelsSoybeans + (1|CountyName),data=landsat[-33,])
  newXandZ <- unique(data.frame(PixelsCorn=landsat$MeanPixelsCorn,
                                PixelsSoybeans=landsat$MeanPixelsSoybeans,
                                CountyName=landsat$CountyName))
  (cpredvar <- get_cPredVar(fitCorn, newdata=newXandZ, nsim=200L, seed=123)) # serial computation
  (cpredvar <- get_cPredVar(fitCorn, newdata=newXandZ, nsim=200L, seed=123, 
        nb_cores=parallel::detectCores()-1L, fit_env=list2env(list(newXandZ=newXandZ))))
}

# GLMM example from Booth & Hobert 1998 JASA
npos <- c(11,16,14,2,6,1,1,4,10,22,7,1,0,0,1,6)
ntot <- c(36,20,19,16,17,11,5,6,37,32,19,17,12,10,9,7)
treatment <- c(rep(1,8),rep(0,8))
clinic <-c(seq(8),seq(8))
clinics <- data.frame(npos=npos,nneg=ntot-npos,treatment=treatment,clinic=clinic)
#
fitClinics <- HLfit(cbind(npos,nneg)~treatment+(1|clinic),family=binomial(),data=clinics)
#
(get_cPredVar(fitClinics, newdata=clinics[1:8,], nsim=200L, seed=123))  # serial computation
(get_cPredVar(fitClinics, newdata=clinics[1:8,], nsim=200L, seed=123, 
      nb_cores=parallel::detectCores()-1, fit_env=list2env(list(clinics=clinics))))
}
}

