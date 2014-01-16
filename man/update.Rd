\name{update.HLfit}
%- \alias{update.HLCor}
\alias{update.HLfit}
\alias{update}
\title{
  Updates an HLCor or HLfit fit 
}
\description{
\code{update} will update and (by default) re-fit a model. It does this mostly by extracting the call stored in the object, updating the call and evaluating that call. 
(however, currently the \code{predictor} argument is processed differently).
}
\usage{
%- \method{update}{HLCor}(object,formula.,...)
\method{update}{HLfit}(object, formula., ..., evaluate = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  A return object from an HLfit call.
}

\item{formula.}{Changes to the formula. Beware of the syntax: see \cite{\link{update.formula}} for details.}

\item{\dots}{
Additional arguments to the call, or arguments with changed values. Use name = NULL to remove the argument name.
}

\item{evaluate}{
  If TRUE, evaluate the new call else return the call.
}
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
wFit <- HLfit(y ~X1*X3+X2*X3+I(X2^2)+(1|batch),family=Gamma(log),
          resid.formula = ~ X3+I(X3^2) ,data=wafers)

# For estimates given by Lee et al., Appl. Stochastic Models Bus. Ind. (2011) 27:  315-328:
# Refit with given beta or/and phi values:
 
betavals <- c(5.55,0.08,-0.14,-0.21,-0.08,-0.09,-0.09)
# reconstruct fitted phi value from predictor for log(phi)
Xphi <- with(wafers,cbind(1,X3,X3^2)) ## design matrix
phifit <- exp(Xphi \%*\% c(-2.90,0.1,0.95))
update(wFit,formula.= . ~ offset(wFit$X \%*\% betavals)+(1|batch),
       ranFix=list(lambda=exp(-3.67),phi=phifit))

## There are subtlety in performing REML fits of constrained models: 
update(wFit,formula.= . ~ offset(wFit$X \%*\% betavals)+(1|batch))
## ... changes the REML correction. Consider instead
update(wFit,formula.= . ~ offset(wFit$X \%*\% betavals)+(1|batch),
       REMLformula=wFit$predictor)
## Alternatively, show original wFit as differences from betavals:  
update(wFit,formula.= . ~ . +offset(wFit$X \%*\% betavals))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
