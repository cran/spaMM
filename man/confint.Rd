\name{confint.HLfit}
\alias{confint.HLfit}
\alias{confint}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Confidence intervals
}
\description{
This computes confidence intervals for given parameter(s), based either on parametric bootstrap or, for \bold{fixed-effect} parameters, on the \code{p_v}-based approximation of the profile likelihood ratio for this parameter. The profiling is over all other fitted parameters, including fixed effects, as well as variances of random effects and spatial correlations if these were fitted. The bootstrap is performed if the \code{parm} argument is a function or a quoted expression or if the \code{boot_args} argument is a list. The profile confidence interval is computed if neither of these conditions is true; in that case \code{parm} must be the name(s) of some \bold{fixed-effect} coefficient.
}
\usage{
\method{confint}{HLfit}(object, parm, level=0.95, verbose=TRUE, 
                          boot_args=NULL, format="default", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of class \code{HLfit}, as returned by the fitting functions in \code{spaMM}.}
  \item{parm}{character vector, integer vector, or function, or a quoted expression. If \bold{character}, the name(s) of parameter(s) to be fitted; if \bold{integer}, their position in the \code{fixef(object)} vector. Valid names are those of this vector. If a \bold{function}, it must return a (vector of) parameter estimate(s) from a fit object. If a \bold{quoted expression}, it must likewise extract parameter estimate(s) from a fit object; this expression must refer to the fitted object as \sQuote{hlfit} (see Examples).}
  \item{level}{The coverage of the interval.}
  \item{verbose}{whether to print the interval or not. As the function returns its more extensive results invisibly, 
  this printing is the only visible output.}
  \item{boot_args}{NULL or a list of arguments passed to functions \code{\link{spaMM_boot}} and \code{\link[boot]{boot.ci}}. It must contain element \code{nsim} (for \code{spaMM_boot}). The \code{type} argument of \code{boot.ci} can only be given as element \code{ci_type}, to avoid conflict with the \code{type} argument of \code{spaMM_boot}.
  }
  \item{format}{Only effective non-default value is \code{"stats"} to return results in the format of the \code{stats::\link[stats]{confint}} result (see Value).}
  \item{\dots}{Additional arguments (maybe not used, but conforming to the generic definition of \code{confint}).}
}
\value{
The format of the value varies, but in all cases distinguished below, one or more tables are included, as a \code{table} attribute, in the format of the \code{stats::\link[stats]{confint}} result, to facilitate consistent extraction of results. By default \code{confint} returns invisibly the full values described below, but if \code{format="stats"}, only the \code{table} attribute is returned. 

If a profile CI has been computed for a single parameter, a list is returned including the confidence \code{interval} as shown by \code{verbose=TRUE}, and the fits \code{lowerfit} and \code{upperfit} giving the profile fits at the confidence bounds. This list bears the \code{table} attribute. 

If a profile CI has been computed for several parameters, a structured list, named according to the parameter names, of such single-parameter results is returned, and a single \code{table} attribute for all parameters is attached to the upper level of the list. 

If a bootstrap was performed, for a single parameter the result of the \code{boot.ci} call is returned, to which a \code{table} attribute is added. This \code{table} is now a list of tables for the different bootstrap CI types (default being \code{normal}, \code{percent}, and \code{basic}), each such table in the format of the \code{stats::confint} results. For several parameters, a named list of \code{boot.ci} results is returned, its names being the parameter names, and the \code{table} attribute is attached to the upper level of the list. 

The \code{boot.ci} return value for each parameter includes the call to \code{boot.ci}. This call is typically shown including a long \code{t} vector, which makes a bulky display. spaMM hacks the printing to abbreviate long \code{t}s.
}
\examples{
% checked in test-confint.R
data("wafers")
wfit <- HLfit(y ~X1+(1|batch), family=Gamma(log), data=wafers, method="ML")
confint(wfit,"X1")  # profile CI
if (spaMM.getOption("example_maxtime")>30) {
   
   # bootstrap CI induced by 'boot_args':
   confint(wfit,names(fixef(wfit)), boot_args=list(nsim=99, seed=123)) 
   
   # bootstrap CI induced by 'parm' being a function:
   confint(wfit,parm=function(v) fixef(v), 
           boot_args=list(nb_cores=10, nsim=199, seed=123))
   
   # Same effect if 'parm' is a quoted expression in terms of 'hlfit':
   confint(wfit,parm=quote(fixef(hlfit)), 
           boot_args=list(nb_cores=10, nsim=199, seed=123))
           
   # CI for the variance of the random effect:          
   ci <- confint(wfit,parm=function(fit){get_ranPars(fit)$lambda[1]}, 
        boot_args=list(nb_cores=10, nsim=199, seed=123))
   # The distribution of bootstrap replicates:
   plot(ecdf(ci$call$t))
   # We may be far from ideal condition for accuracy of bootstrap intervals;
   # for variances, a log transformation may sometimes help, but not here.
   
   # Passing arguments to child processes, as documented in help("spaMM_boot"):
   set.seed(123)
   rvar <- runif(nrow(wafers))
   wrfit <- fitme(y ~X1+(1|batch), family=Gamma(log), data=wafers, fixed=list(phi=rvar))
   confint(wrfit, parm = "(Intercept)", boot_args = list(nsim = 100, nb_cores = 2,
           fit_env = list(rvar=rvar)))
   
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
