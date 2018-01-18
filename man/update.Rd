\name{update.HLfit}
%- \alias{update.HLCor}
\alias{update.HLfit}
\alias{update}
\title{
  Updates an HLCor or HLfit fit 
}
\description{
\code{update} will update and (by default) re-fit a model. It does this mostly by extracting the call stored in the object, updating the call and evaluating that call. 
(however, currently the \code{predictor} argument is processed differently). Using \code{update} is a risky programming style (see Note). 
}
\note{
\code{update}, as a general rule, is tricky. For example with \code{m1 <- glm(y~x + offset(<...>),<...>); m2 <- update(m1,.~.-x)}, \code{m2} will lose the intercept term, which would be retained in the absence of an offset. \code{update.HLfit} has the same behaviour. Further, \code{update} is easily affected in a non-transparent way by changes in variables used in the original call. For example
\code{
foo <- rep(1,10)
m <- lm(rnorm(10)~1, weights=foo)
rm(foo)
update(m, .~.) # Error
}
To avoid such problems, \code{spaMM} tries to avoid references to variables in the global environment, by enforcing that the data are explicitly provided by the \code{data} argument, and that prior weights, if not constant, are in the data.
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
data("wafers")
## First the fit to be updated:
wFit <- HLfit(y ~X1*X3+X2*X3+I(X2^2)+(1|batch),family=Gamma(log),
          resid.model = ~ X3+I(X3^2) ,data=wafers)

# For estimates given by Lee et al., Appl. Stochastic Models Bus. Ind. (2011) 27:  315-328:
# Refit with given beta or/and phi values:
 
betavals <- c(5.55,0.08,-0.14,-0.21,-0.08,-0.09,-0.09)
# reconstruct fitted phi value from predictor for log(phi)
Xphi <- with(wafers,cbind(1,X3,X3^2)) ## design matrix
phifit <- exp(Xphi \%*\% c(-2.90,0.1,0.95))
update(wFit,formula.= . ~ offset(wFit$`X.pv` \%*\% betavals)+(1|batch),
       ranFix=list(lambda=exp(-3.67),phi=phifit))

## There are subtlety in performing REML fits of constrained models,
##   illustrated by thefact thatthe following fir does not recover 
##   the original likelihood values, because dispersion parameters are
##   estimated but the REML correction changes with the formula:
update(wFit,formula.= . ~ offset(wFit$`X.pv` \%*\% fixef(wFit))+(1|batch))
## To maintain the original REML correction, Consider instead
update(wFit,formula.= . ~ offset(wFit$`X.pv` \%*\% fixef(wFit))+(1|batch),
       REMLformula=wFit$predictor)  ## recover original p_v and p_bv 
## Alternatively, show original wFit as differences from betavals:  
update(wFit,formula.= . ~ . +offset(wFit$`X.pv` \%*\% betavals))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
