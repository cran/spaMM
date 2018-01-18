\name{fixed}
\alias{fixed}
\alias{ranFix}
\alias{ranPars}
\alias{etaFix}
\alias{ranCoefs}

\title{Fixing some parameters}
\description{The fitting functions allow some parameters to be fixed rather than estimated, by way of the \code{fixed} argument for \code{fitme}, of \code{etaFix} (linear predictor coefficents) and \code{ranFix} (random-effect and residual dispersion parameters) for \code{HLfit} and \code{corrHLfit}, and of \code{ranPars} for \code{HLCor}. The diversity of names may be confusing, but keep in mind that \code{ranFix} allows one to fix parameters that \code{HLfit} and \code{corrHLfit} would otherwise estimate, while \code{ranPars} can be used to set required parameters for \code{HLCor}, which it would otherwise be unable to estimate (e.g., Matern correlation parameters). 

Each of these arguments is a list.

\code{ranFix} elements taken into account by \code{HLfit} include \code{phi} (variance of residual error, for gaussian and Gamma HGLMs), \code{lambda} (random-effect variances), and \code{ranCoefs} (variance-correlation information for random-coefficient terms). To assign values for only some random-effect terms, \code{lambda} and \code{ranCoefs} can be incomplete, e.g. \code{lambda=c(NA,1)} or \code{lambda=c("2"=1)} (note the name) to assign a value only to the variance of the second of two random effects.
\code{ranCoefs} is a \code{list} of numeric vectors, each numeric vector specifying the variance and correlation parameters for a random-coefficient term. This input matches the printed summary of a fit. The elements must be given in the order of the \code{\link{lower.tri}} of a covariance matrix, as shown e.g. by\cr
  \code{m2 <- matrix(NA, ncol=2,nrow=2); m2[lower.tri(m2,diag=TRUE)] <- seq(3); m2}.\cr 
  For example, to assign variances values 3 and 7, and correlation value -0.05, to a second random effect, one can use \code{ranCoefs=list("2"=c(3,-0.05,7))} (note the name).
  
Additional \code{ranFix} elements are taken into account by \code{corrHLfit}, as follows.
  For the Matérn model, these are the correlation parameters are \code{rho} (scale parameter(s)), \code{nu} (smoothness parameter), and (optional) \code{Nugget} (see \code{\link{Matern}}). The \code{rho} parameter can itself be a vector with different values for different geographic coordinates.
  For the \code{adjacency} model, the only correlation parameter is a scalar \code{rho} (see \code{\link{adjacency}}).
  For the \code{AR1} model, the only correlation parameter is a scalar \code{ARphi} (see \code{\link{AR1}}).

\code{ranPars} elements taken into account by \code{HLCor} include all the above \code{ranFix} elements.

The only \code{etaFix} element considered here is \code{beta}, which should be a vector of (a subset of) the coefficients (\eqn{\beta}) of the fixed effects, with names as shown in a fit without such given values. In contrast to an offset specification, it affects by default the REML correction for estimation of dispersion parameters, which depends only on which \eqn{\beta} coefficients are estimated rather than given. This default behaviour will be overridden whenever a non-null \code{REMLformula} is provided to \code{HLfit} or the other fitting functions (see Example).

\code{fixed} elements taken into account by \code{fitme} include all the above \code{ranFix} and \code{etaFix} elements.

}
\details{
\code{REMLformula} is the preferred way to control non-standard REML fits. Alternatively, with a non-NULL \code{etaFix$beta}, REML can also be performed as if all \eqn{\beta} coefficients were estimated, by adding attribute \code{keepInREML=TRUE} to \code{etaFix$beta}. Using an  \code{REMLformula} will override such a specification.  
}
\examples{
\dontrun{
data("wafers")
# Fixing random-coefficient parameters:
HLfit(y~X1+(X2|batch), data=wafers, ranFix=list(ranCoefs=list("1"=c(2760, -0.1, 1844))))
# fixing coefficients of the linear predictor, but with REML as if they were not fixed:
HLfit(y ~X1+X2+X1*X3+X2*X3+I(X2^2)+(1|batch), data=wafers, family=Gamma(log),
      etaFix=list(beta=c("(Intercept)"=5.61208, X1=0.08818, X2=-0.21163, X3=-0.13948, 
                         "I(X2^2)"=-0.10378, "X1:X3"=-0.08987, "X2:X3"=-0.08779)),
      REMLformula=y ~X1+X2+X1*X3+X2*X3+I(X2^2)+(1|batch))

data("Loaloa")
# Fixing some Matern correlation parameters, in corrHLfit:
corrHLfit(cbind(npos,ntot-npos)~elev1+elev2+elev3+elev4+maxNDVI1+seNDVI
           +Matern(1|longitude+latitude),
             data=Loaloa,family=binomial(),ranFix=list(nu=0.5,Nugget=2/7))
# Fixing all mandatory Matern correlation parameters, in HLCor:
HLCor(cbind(npos,ntot-npos)~elev1+elev2+elev3+elev4+maxNDVI1+seNDVI
           +Matern(1|longitude+latitude),
             data=Loaloa,family=binomial(),ranPars=list(nu=0.5,rho=0.7))
}
}

