\name{update.HLfit}
\alias{update.HLfit}
\alias{update}
\alias{update_resp}
\alias{refit}
\title{
  Updates a fit 
}
\description{
\code{update} and \code{update_resp} will update and (by default) re-fit a model. They do this mostly by extracting the call stored in the object, updating the call and evaluating that call. Using \code{update} is a risky programming style (see Note). 
\code{update_resp} handles a new response vector as produced by \code{simulate}. 
}
\note{
\code{update}, as a general rule, is tricky. \code{update} methods are easily affected in a non-transparent way by changes in variables used in the original call. For example
\code{
foo <- rep(1,10)
m <- lm(rnorm(10)~1, weights=foo)
rm(foo)
update(m, .~.) # Error
}
To avoid such problems, \code{spaMM} tries to avoid references to variables in the global environment, by enforcing that the data are explicitly provided to the fitting functions by the \code{data} argument, and that any variable used in the \code{prior.weights} argument is in the data.

\code{spaMM}'s \code{update} method was all the more tricky when spaMM called \code{stats::update.formula} whose results endorse \code{stats}'s (sometimes annoying) convention that a formula without an explicit intercept term actually includes an intercept. \code{spaMM::update} now avoids this problem. \bold{Formula updates should still be carefully checked}, as getting them perfect has not been on the priority list.
}
\usage{
\method{update}{HLfit}(object, formula., ..., evaluate = TRUE)
update_resp(object, newresp, ..., evaluate = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  A return object from an HLfit call.
}

\item{formula.}{Changes to the formula. Beware of the syntax: see \cite{\link{update.formula}} for details.}
\item{newresp}{New response vector.}

\item{\dots}{
Additional arguments to the call, or arguments with changed values. Use \emph{name}\code{ = NULL} to remove the argument with given \emph{name}.
}

\item{evaluate}{
  If TRUE, evaluate the new call else return the call.
}
}
\value{
  An HLfit fit of the same type as the input object, or a call object, depending on the \code{evaluate} value.
}
\seealso{
See also \code{\link{HLCor}},  \code{\link{HLfit}}.
}
\examples{
data("wafers")
## First the fit to be updated:
wFit <- HLfit(y ~X1*X3+X2*X3+I(X2^2)+(1|batch),family=Gamma(log),
          resid.model = ~ X3+I(X3^2) ,data=wafers)

newresp <- simulate(wFit)
update_resp(wFit,newresp=newresp)

# For estimates given by Lee et al., Appl. Stochastic Models Bus. Ind. (2011) 27:  315-328:
# Refit with given beta or/and phi values:
 
betavals <- c(5.55,0.08,-0.14,-0.21,-0.08,-0.09,-0.09)
# reconstruct fitted phi value from predictor for log(phi)
Xphi <- with(wafers,cbind(1,X3,X3^2)) ## design matrix
phifit <- exp(Xphi \%*\% c(-2.90,0.1,0.95))
upd_wafers <- wafers
upd_wafers$off_b <- wFit$`X.pv` \%*\% betavals
update(wFit,formula.= . ~ offset(off_b)+(1|batch), data=upd_wafers,
       ranFix=list(lambda=exp(-3.67),phi=phifit))

## There are subtlety in performing REML fits of constrained models,
##   illustrated by the fact that the following fit does not recover 
##   the original likelihood values, because dispersion parameters are
##   estimated but the REML correction changes with the formula:
upd_wafers$off_f <- wFit$`X.pv` \%*\% fixef(wFit) ## = predict(wFit,re.form=NA,type="link")
update(wFit,formula.= . ~ offset(off_f)+(1|batch), data=upd_wafers)
## To maintain the original REML correction, Consider instead
update(wFit,formula.= . ~ offset(off_f)+(1|batch), data=upd_wafers,
       REMLformula=formula(wFit))  ## recover original p_v and p_bv     
## Alternatively, show original wFit as differences from betavals:  
update(wFit,formula.= . ~ . +offset(off_f), data=upd_wafers)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
