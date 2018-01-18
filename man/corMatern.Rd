\name{corMatern}
\alias{corMatern}
\alias{Initialize.corMatern} 
\alias{Variogram.corMatern} 
\alias{coef.corMatern} 
\alias{coef<-.corMatern} 
\alias{corFactor.corMatern} 
\alias{corMatrix.corMatern} 
\alias{getCovariate.corMatern} 
\alias{logDet.corMatern} 
\alias{recalc.corMatern} 
\title{Matern Correlation Structure as a corSpatial object}
\description{
  This implements the Matérn correlation structure (see \code{\link{Matern}}) for use with \code{lme} or \code{glmmPQL}. Usage is as for others \code{corSpatial} objects 
  such as \code{corGaus} or \code{corExp}, except that the Matérn family
  has an additional parameter. This function was defined for comparing
  results obtained with \code{corrHLfit} to those produced by \code{lme} and \code{glmmmPQL}. There
  are problems in fitting (G)LMMs in the latter way, so it is not a
  recommended practice. 
}
\usage{
corMatern(value = c(1, 0.5), form = ~1, nugget = FALSE, nuScaled = FALSE, 
          metric = c("euclidean", "maximum", "manhattan"), fixed = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{value}{
   An optional vector of parameter values, with serves as initial values or as fixed values depending on the \code{fixed} argument. 
   It has either two or three elements, depending on the \code{nugget} argument.

If \code{nugget} is \code{FALSE}, \code{value} should have
   two elements, corresponding to the "range" and the "smoothness" \eqn{\nu} of the
   Matérn correlation structure. If \code{value} has zero length, 
   the default is a range of 90\% of the minimum distance and a
   smoothness of 0.5 (exponential correlation). 
   \bold{Warning}: the range parameter used in \code{corSpatial} objects is the inverse of the scale 
   parameter used in \code{\link{MaternCorr}} and thus they have opposite meaning despite both being denoted \eqn{\rho} elsewhere in this package or in \code{nlme} literature. 
   
   If \code{nugget} is \code{TRUE}, meaning that a nugget effect
   is present, \code{value} can contain two or three elements, the first
   two as above, the third  being the "nugget effect" (one minus the
   correlation between two observations taken arbitrarily close
   together). If \code{value} has length zero or two, the nugget defaults to 0.1. 
   The range and smoothness must be greater than zero and the nugget must be
   between zero and one. 
  }
  \item{form}{
   (Pasted from corSpatial) a one sided formula of the form \code{~ S1+...+Sp}, or
   \code{~ S1+...+Sp | g}, specifying spatial covariates \code{S1}
   through \code{Sp} and,  optionally, a grouping factor \code{g}. 
   When a grouping factor is present in \code{form}, the correlation
   structure is assumed to apply only to observations within the same
   grouping level; observations with different grouping levels are
   assumed to be uncorrelated. Defaults to \code{~ 1}, which corresponds
   to using the order of the observations in the data as a covariate,
   and no groups.}
  \item{nugget}{
   an optional logical value indicating whether a nugget
   effect is present. Defaults to \code{FALSE}.
}
  \item{nuScaled}{
    If \code{nuScaled} is set to \code{TRUE} the "range" parameter
    \eqn{\rho} is divided by \eqn{2 \sqrt\nu}. With this option and
    for large values of \eqn{\nu}, \code{corMatern} reproduces the
    calculation of \code{corGaus}. Defaults to \code{FALSE}, in which
    case the function compares to \code{corGaus} with range parameter
    \eqn{2(\sqrt\nu)\rho} when \eqn{\nu} is large.
}
  \item{metric}{
   (Pasted from corSpatial) an optional character string specifying the distance
   metric to be used. The currently available options are
   \code{"euclidean"} for the root sum-of-squares of distances;
   \code{"maximum"} for the maximum difference; and \code{"manhattan"}
   for the sum of the absolute differences. Partial matching of
   arguments is used, so only the first three characters need to be
   provided. Defaults to \code{"euclidean"}.
}
  \item{fixed}{
   an optional logical value indicating whether the
   coefficients should be allowed to vary in the optimization, or kept
   fixed at their initial value. Defaults to \code{FALSE}, in which case
   the coefficients are allowed to vary.
}
}
\details{
  This function is a constructor for the \code{corMatern} class,
  representing a Matérn spatial correlation structure. See
   \code{\link{MaternCorr}} for details on the Matérn family. 
}
\note{The R and C code for the methods for \code{corMatern} objects builds on code for \code{corSpatial} objects, by D.M. Bates, J.C. Pinheiro and S. DebRoy, in a circa-2012 version of nlme.}
\value{
  an object of class \code{corMatern}, also inheriting from class
  \code{corSpatial}, representing a Matérn spatial correlation
  structure.
}
\references{
Mixed-Effects Models in S and S-PLUS, José C. Pinheiro and Douglas M. Bates, Statistics and Computing Series, Springer-Verlag, New York, NY, 2000.
}
\seealso{
  \code{\link{glmmPQL}}, \code{\link{lme}}
}
\examples{
## LMM
data("blackcap")
blackcapD <- cbind(blackcap,dummy=1) ## obscure, isn't it? 
## With method= 'ML' in lme, The correlated random effect is described 
##  as a correlated residual error and no extra residual variance is fitted:
nlme::lme(fixed = migStatus ~ means, data = blackcapD, random = ~ 1 | dummy, 
    correlation = corMatern(form = ~ longitude+latitude  | dummy), 
    method = "ML")

## Binomial GLMM
if (spaMM.getOption("example_maxtime")>32) {
 data("Loaloa")
 LoaloaD <- cbind(Loaloa,dummy=1) 
 MASS::glmmPQL(fixed =cbind(npos,ntot-npos)~elev1+elev2+elev3+elev4+maxNDVI1+seNDVI, 
        data = LoaloaD, random = ~ 1 | dummy,family=binomial, 
        correlation = corMatern(form = ~ longitude+latitude | dummy))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
