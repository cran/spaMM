\name{confint.HLfit}
\alias{confint.HLfit}
\alias{confint}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Confidence intervals
}
\description{
This computes confidence intervals for a given parameter, based either on parametric bootstrap or, for \bold{fixed-effect} parameters, on the \code{p_v}-based approximation of the profile likelihood ratio for this parameter. The profiling is over all other fitted parameters: other fixed effects, as well as variances of random effects and spatial correlations if these were fitted. The bootstrap is performed if the \code{parm} argument is a function or a quoted expression or if the \code{boot_args} argument is a list. The profile confidence interval is computed if neither of these conditions is true; in that case \code{parm} must be the name of a \bold{fixed-effect} coefficient.
}
\usage{
\method{confint}{HLfit}(object, parm, level=0.95, verbose=TRUE, 
                          boot_args=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of class \code{HLfit}, as returned by the fitting functions in \code{spaMM}.}
  \item{parm}{character, or function, or a quoted expression. If \bold{character}, The name of a parameter to be fitted, or its position in the the object's \code{$fixef} vector. Valid names are those of the object's \code{$fixef}. If a \bold{function}, it must return the focal parameter estimate from a fit object. If a \bold{quoted expression}, it must likewise extract a parameter estimate from a fit object; this expression must refer to the fitted object as \sQuote{hlfit} (see Examples).}
  \item{level}{The coverage of the interval.}
  \item{verbose}{whether to print the interval or not. As the function returns its more extensive results invisibly, 
  this printing is the only visible output.}
  \item{boot_args}{NULL or a list of arguments passed to functions \code{\link{spaMM_boot}} and \code{\link[boot]{boot.ci}}. It must contain element \code{nsim} (for \code{spaMM_boot}). The \code{type} argument of \code{boot.ci} can only be given as element \code{ci_type}, to avoid conflict with the \code{type} argument of \code{spaMM_boot}.
  }
  \item{\dots}{Additional arguments (maybe not used, but conforming to the generic definition of \code{confint}).}
}
\value{If a bootstrap was performed, the result of the \code{boot.ci} call is returned, with basic and percentile intervals, and modified \code{call} element. Otherwise, a list including the confidence \code{interval} for the target parameter, and the fits \code{lowerfit} and \code{upperfit} giving the profile fits at the confidence bounds. This is returned invisibly.
}
\examples{
% checked in test-confint.R
data("wafers")
wfit <- HLfit(y ~X1+(1|batch),family=Gamma(log),data=wafers,method="ML")
confint(wfit,"X1")  # profile CI
if (spaMM.getOption("example_maxtime")>30) {
   confint(wfit,"X1", boot_args=list(nsim=99)) # bootstrap CI induced by 'boot_args'
   # bootstrap CI induced by 'parm' being a function:
   confint(wfit,parm=function(v) fixef(v)[["X1"]], 
           boot_args=list(nb_cores=10, nsim=199))
   # Same effect if 'parm' is a quoted expression in terms of 'hlfit':
   confint(wfit,parm=quote(fixef(hlfit)[["X1"]]), 
           boot_args=list(nb_cores=10, nsim=199))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
