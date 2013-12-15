\name{update.HL}
\alias{update.HL}
\alias{update}
\alias{update.HLCor}
\alias{update.HLfit}
\title{
  Updates an HLCor or HLfit fit 
}
\description{
  This is a convenient interface for updating a fit by only specifiying the arguments that are modified in the new fit.
}
\usage{
\method{update}{HLCor}(object, ...)
\method{update}{HLfit}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  A return object from an HLCor or HLfit call.
}
  \item{\dots}{
  New HLCor/HLfit arguments
}
}
\details{
  In principle any HLCor or HLfit argument could be updated and the model refit. This function is useful in particular to refit a model for given parameter values. 
  However, currently this is straightforward only for complete parameter vectors, e.g. a complete \eqn{\beta} vector which can be provided as an element of \code{etaFix}
 (see Examples). 
  
  To specify values for an incomplete \eqn{\beta}, one should not use the \code{etaFix} argument. Rather, one needs to provide the resulting offset values 
  and change the model formula to exclude the corresponding predictor variables (which is what is done automatically by \code{update} in the case of an \code{etaFix} update).

  The same should apply to \eqn{\phi} in the case a model formula has been given for \eqn{\phi}. 
  It is also possible to provide the \eqn{\phi} values as \code{ranFix$phi} (see Examples). 
}
\value{
  An HLCor or HLfit fit of the same type as the input object.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also \code{\link{HLCor}},  \code{\link{HLfit}}.
}
\examples{
data(wafers)
## First the fit to be updated:
wFit <- HLfit(predictor= y ~X1*X3+X2*X3+I(X2^2)+(1|batch),family=Gamma(log),
          resid.predictor = ~ X3+I(X3^2) ,data=wafers)

# For estimates given by Lee et al., Appl. Stochastic Models Bus. Ind. (2011) 27:  315-328:
# Refit with given beta or/and phi values:
 
betavals <- c(5.55,0.08,-0.14,-0.21,-0.08,-0.09,-0.09)
update(wFit,etaFix=list(beta=betavals))
# (p_bv is not comparable to the previous p_bv as the fixed-effects model is different).

# reconstruct fitted phi value from predictor for log(phi)
Xphi <- with(wafers,cbind(1,X3,X3^2)) ## design matrix
phifit <- exp(Xphi \%*\% c(-2.90,0.1,0.95))
update(wFit,ranFix=list(phi=exp(Xphi \%*\% c(-2.90,0.1,0.95))))
update(wFit,etaFix=list(beta=betavals),ranFix=list(lambda=exp(-3.67),phi=phifit))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
