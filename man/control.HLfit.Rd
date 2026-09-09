\name{control.HLfit}
\alias{control.HLfit}
\title{
Control parameters of the HLfit fitting algorithm 
}
\description{
  A list of parameters controlling the \code{\link{HLfit}} fitting algorithm (called by all fitting functions in \pkg{spaMM}). It can mostly be ignored in routine use.
}
\usage{
# <fitting function>(., control.HLfit=list(...)) 
}
\details{
  Possible controls are:
  
\itemize{
  \item{\code{algebra, sparse_precision}}: see \code{\link{algebra}}.
  
  \item{\code{max.iter}}: the number of iterations of the iterative algorithm for joint estimation of dispersion parameters and of coefficients of the linear predictor. Defaults to 200. This is typically much more than necessary, unless there is little information to separately estimate \eqn{\lambda} and \eqn{\phi} parameters.
  
  \item{ \code{conv.threshold} and \code{spaMM_tol}}: \code{spaMM_tol} is a list of tolerance values, with elements \code{Xtol_rel} and \code{Xtol_abs} that define thresholds for relative and absolute changes in parameter values in iterative algorithms (used in tests of the form \dQuote{d(param)< Xtol_rel * param + Xtol_abs}, so that \code{Xtol_abs} is operative only for small parameter values). \code{conv.threshold} is the older way to control \code{Xtol_rel}. Default values are given by \code{spaMM.getOption("spaMM_tol")}.
  
  \item{ \code{min_phi}}: NULL or a positive numeric value, being the minimum value of the residual dispersion parameter \eqn{\phi}, when \dQuote{inner-estimated} (\code{\link{inner-vs-outer}} for this concept). Values approaching zero may induce numerical errors in various steps of the fitting algorithm. For Gaussian- or Gamma-family models, \eqn{\phi} values as low as \code{1e-12} are usually not a source of concern; but setting \code{min_phi} to higher values has been needed at least for fits with an experimental implementation of the \code{\link{tweedie}} family.
  
  \item{ \code{break_conv_logL}}: a boolean specifying whether the iterative algorithm should terminate when log-likelihood appears to have converged (roughly, when its relative variation over on iteration is lower than 1e-8). Default is FALSE (convergence is then assessed on the parameter estimates rather than on log-likelihood).
  
  \item{\code{iter.mean.dispFix}}: the number of iterations of the iterative algorithm for coefficients of the linear predictor, if no dispersion parameters are estimated by the iterative algorithm. Defaults to 200 except for Gamma(log)-family models.
  
  \item{\code{iter.mean.dispVar}}: the number of iterations of the iterative algorithm for coefficients of the linear predictor, if some dispersion parameter(s) is estimated by the iterative algorithm. Defaults to 50 except for Gamma(log)-family models.  
}

\code{resid.family} was a previously documented control (before version 2.6.40), and might still operate as previously documented, but should not be used in new code.
  %%%%%  As it just a list member and not formal arg, we can just remove the doc.
%  \code{resid.family} Allows one to change the link for modeling of residual variance \eqn{\phi}, which is \code{"log"} by default. The family is always Gamma, so the non-default possible values of \code{resid.family} are \code{Gamma(identity)} or \code{Gamma(inverse)}. Only the default value ensures that the fitted \eqn{\phi} is positive.     

}
