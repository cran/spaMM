\name{fixed}
\alias{fixed}
\alias{ranFix}
\alias{ranPars}
\alias{etaFix}
\alias{ranCoefs}
\alias{corrPars}

\title{Fixing some parameters}
\description{The fitting functions allow some parameters to be fixed rather than estimated. Fixed-effect coefficients can be set by by way of the \code{etaFix} argument (linear predictor coefficents) for all fitting functions. Random-effect parameters can be set by an argument with a different name for the different fitting functions: \code{fixed} for \code{fitme}, \code{ranFix} for \code{HLfit} and \code{corrHLfit}, and \code{ranPars} for \code{HLCor}. This diversity of names may be confusing, but keep in mind that \code{ranFix} allows one to fix parameters that \code{HLfit} and \code{corrHLfit} would otherwise estimate, while \code{ranPars} can be used to set correlation parameters that \code{HLCor} does not estimate but nevertheless requires (e.g., Matérn parameters).

}
\details{
\bold{etaFix} is a list with principal element \code{beta}, which should be a vector of (a subset of) the coefficients (\eqn{\beta}) of the fixed effects, with names as shown in a fit without such given values. If REML is used to fit random effect parameters, then \code{etaFix} affects by default the REML correction for estimation of dispersion parameters, which depends only on which \eqn{\beta} coefficients are estimated rather than given. This default behaviour will be overridden whenever a non-null \code{REMLformula} is provided to the fitting functions (see Example). \code{REMLformula} is the preferred way to control non-standard REML fits. Alternatively, with a non-NULL \code{etaFix$beta}, REML can also be performed as if all \eqn{\beta} coefficients were estimated, by adding attribute \code{keepInREML=TRUE} to \code{etaFix$beta}. Using an  \code{REMLformula} will override such a specification.  

Despite its different name for different fitting functions, the argument for fixing random-effect parameters has a common syntax for all functions. It is a list, with the following possible elements, whose nature is further detailed below: \bold{phi} (variance of residual error, for gaussian and Gamma HGLMs), \bold{lambda} (random-effect variances, except for random-coefficient terms), \bold{ranCoefs} (random-coefficient parameters), and \bold{corrPars} (correlation parameters, when handled by the fitting function). Individual correlation parameters such as \bold{rho, nu, Nugget, ARphi}... are also possible When there is no ambiguity as to which random effect these correlation parameters apply. This was the original syntax, conceived when \code{spaMM} handled a single spatial random effect, and it is still convenient when applicable, but it should not be mixed with \code{corrPars} usage. 

\bold{phi} may be a single value or a vector of the same length as the response vector (the number of rows in the \code{data}, once non-informative rows are removed). 

\bold{lambda} may be a single value (if there is a single random effect, or a vector allowing to specify unambiguously variance  values for some random effect(s). It can thus take the form \code{lambda=c(NA,1)} or \code{lambda=c("2"=1)} (note the name) to assign a value only to the variance of the second of two random effects.

\bold{ranCoefs} is a \code{list} of numeric vectors, each numeric vector specifying the variance and correlation parameters for a random-coefficient term. As for \code{lambda}, it may be incomplete, using names to specify the random effect to which the parameters apply. For example, to assign variances values 3 and 7, and correlation value -0.05, to a second random effect, one can use \code{ranCoefs=list("2"=c(3,-0.05,7))} (note the name). The elements of each vector are variances and correlations, matching those of the printed summary of a fit. The order of these elements must be the order of the \code{\link{lower.tri}} of a covariance matrix, as shown e.g. by\cr
  \code{m2 <- matrix(NA, ncol=2,nrow=2); m2[lower.tri(m2,diag=TRUE)] <- seq(3); m2}.\cr 
\code{fitme} accepts partially fixed parameters for a random coefficient term, e.g., \code{ranCoefs=list("2"=c(NA,-0.05,NA))}, although this may not mix well with some obscure options, such as \code{control=list(refit=list(ranCoefs=TRUE))} which will ignore the fixed values.    
  
\bold{corrPars} is a list, and it may also be incomplete, using names to specify the affected random effect as shown for \code{lambda} and \code{ranCoefs}. For example, \code{ranFix=list(corrPars=list("1"=list(nu=0.5)))} makes explicit that \code{nu=0.5} applies to the first (\code{"1"}) random effect in the model formula. Its elements may be the correlation parameters of the given random effect. For the Matérn model, these are the correlation parameters \code{rho} (scale parameter(s)), \code{nu} (smoothness parameter), and (optionally) \code{Nugget} (see \code{\link{Matern}}). The \code{rho} parameter can itself be a vector with different values for different geographic coordinates. %
  For the \code{adjacency} model, the only correlation parameter is a scalar \code{rho} (see \code{\link{adjacency}}).
  For the \code{AR1} model, the only correlation parameter is a scalar \code{ARphi} (see \code{\link{AR1}}).
  Consult the documentation for other types of random effects, such as \code{\link{Cauchy}} or \code{\link{IMRF}}, for any information missing here. 

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

