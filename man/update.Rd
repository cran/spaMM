\name{update.HLfit}
\alias{update.HLfit}
\alias{update.formula}
\alias{update_formulas}
\alias{update_resp}
\alias{respName}
\alias{refit} % a refit() fn could conflict with the generic defined by lme4
\title{
  Updates a fit 
}
\description{
\code{update} and \code{update_resp} will update and (by default) re-fit a model. They do this mostly by extracting the call stored in the object, updating the call and evaluating that call. Using \code{update(<fit>)} is a risky programming style (see Details). \code{update_formulas(<mv fit>, ...)} can update formulas from a \code{fitmv} fit as well as the single formula of a fit by the other fitting functions.

\code{update_resp} handles a new response vector as produced by \code{simulate}. 
}
\details{
%\code{update_resp} updates various components of the fitted object, but not its \code{data} element, as the latter approach would meet difficulties when the response is an expression (such as \code{I(<some variable>^2)}) rather than simply a variable in the data. This implies that the object returned by \code{update_resp} cannot be used as if the \code{data} themselves had been updated. % _F I X M E__ Adding a protection, e.g. NA's in the data could be nice, but is not straightforward, for he same reason.
\bold{Controlling response updating}: Updating the data may be tricky when the response specified in a formula is not simply the name of a variable in the data. For example,
if the response was specified as \code{I(foo^2)} the variable \code{foo} is not what \code{simulate.HLfit} will simulate, so \code{foo} should not be updated with such simulation results, yet this is what should be updated in the \code{data}. For some time \pkg{spaMM} has handled such cases by using an alternative way to provide updated response information, but this has some limitations. So \pkg{spaMM} now update the data after checking that this is correct, which the consequence that when response updating is needed (notably, for bootstrap procedures), the response should preferably be specified as the name of a variable in the data, rather than a more complicated expression. 

However, in some cases, dynamic evaluation of the response variable may be helpful. For example, for bootstrapping hurdle models, the zero-truncated response may be specified as \code{I(count[presence>0] <- NA; count)} (where both the zero-truncated \code{count} and binary \code{presence} variables are both updated by the bootstrap simulation). In that case the names of the two variables to be updated is provided by setting (say)\cr
\code{<fit object>$respNames <- c("presence", "count")}\cr
for an hurdle model fit as a bivariate-response model, with first submodel for presence/absence, and second submodel for zero-truncated response. A full example is developed in the \dQuote{Gentle introduction} to \pkg{spaMM} (
\url{https://gitlab.mbb.univ-montp2.fr/francois/spamm-ref/-/blob/master/vignettePlus/spaMMintro.pdf}).
Alternatively for univariate-response fits, use \cr 
\code{<fit object>$respName <- "count"}

\bold{Controlling formula updating}: Early versions of \pkg{spaMM}'s \code{update} method relied on \code{stats::update.formula} whose results endorse \code{stats}'s (sometimes annoying) convention that a formula without an explicit intercept term actually includes an intercept. \code{spaMM::update.HLfit} was then defined to avoid this problem. \bold{Formula updates should still be carefully checked}, as getting them perfect has not been on the priority list. 

Various post-fit functions from base R may use \code{update.formula} directly, rather than using automatic method selection for \code{update}. \code{update.formula} is not itself a generic, which leads to the following problem. To make \code{update.formula()} work on multivariate-response fits, one would like to be able to redefine it as a generic, with an \code{HLfit} method that would perform what \code{update_formulas} does, but such a redefinition appears to be forbidden in a package distributed on CRAN. Instead it is suggested to define a new generic \code{spaMM::update}, which could have a \code{spaMM::update.formula} as a method (possibly itself a generic). This would be of limited interest as the new \code{spaMM::update.formula} would be visible to \code{spaMM::update} but not to \code{stats::update}, and thus the post-fit functions from base R would still not use this method.  

\bold{Safe updating}: \code{update(<fit>, ...)}, as a general rule, is tricky. \code{update} methods are easily affected in a non-transparent way by changes in variables used in the original call. For example
\preformatted{
foo <- rep(1,10)
m <- lm(rnorm(10)~1, weights=foo)
rm(foo)
update(m, .~.) # Error
}
To avoid such problems, \pkg{spaMM} tries to avoid references to variables in the global environment, by enforcing that the data are explicitly provided to the fitting functions by the \code{data} argument, and that any variable used in the \code{prior.weights} argument is in the data.

Bugs can also result when calling \code{update} on a fit produced within some function, say function \code{somefn} calling \code{fitme(data=mydata,...)}, as e.g. \code{update(<fit>)} will then seek a global variable \code{mydata} that may differ from the fitted \code{mydata} which was local to \code{somefn}.   
}
\usage{
\method{update}{HLfit}(object, formula., ..., evaluate = TRUE)
update_resp(object, newresp, ..., evaluate = TRUE)

update_formulas(object, formula., ...)

# <fit object>$respName[s] : see Details.
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{object}{
  A return object from an HLfit call.
}
\item{formula.}{
  A standard \code{formula}; or a \code{formula} with a peculiar syntax only describing changes to the original model formula 
  (see \code{\link[stats]{update.formula}} for details); or (for multivariate-response models) a list of \code{formula} of such types.
}
\item{newresp}{New response vector.}
\item{\dots}{
Additional arguments to the call, or arguments with changed values. Use \emph{name}\code{ = NULL} to remove the argument with given \emph{name}.
}

\item{evaluate}{
  If TRUE, evaluate the new call else return the call.
}
}
\value{
  \code{update.formula(object)} returns an object of the same nature as \code{formula(object)}. The other functions and methods return an HLfit fit of the same type as the input object, or a call object, depending on the \code{evaluate} value. \bold{Warning:} The object returned by \code{update_resp} cannot be used safely for further programming, for the reason explained in the Details section.  
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
designX <- get_matrix(wFit)
upd_wafers$off_b <- designX \%*\% betavals
update(wFit,formula.= . ~ offset(off_b)+(1|batch), data=upd_wafers,
       ranFix=list(lambda=exp(-3.67),phi=phifit))

## There are subtlety in performing REML fits of constrained models,
##   illustrated by the fact that the following fit does not recover 
##   the original likelihood values, because dispersion parameters are
##   estimated but the REML correction changes with the formula:
upd_wafers$off_f <- designX \%*\% fixef(wFit) ## = predict(wFit,re.form=NA,type="link")
update(wFit,formula.= . ~ offset(off_f)+(1|batch), data=upd_wafers)
#
## To maintain the original REML correction, Consider instead
update(wFit,formula.= . ~ offset(off_f)+(1|batch), data=upd_wafers,
       REMLformula=formula(wFit))  ## recover original p_v and p_bv     
## Alternatively, show original wFit as differences from betavals:  
update(wFit,formula.= . ~ . +offset(off_f), data=upd_wafers)
}

