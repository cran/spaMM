\name{verbose}
\alias{verbose}

\title{Details of the \code{verbose} argument of various functions}

\description{This partially documents the usage of the \code{verbose} argument in tow contexts: (1) as an argument of the fitting functions; (2) as an  argument of \code{\link{confint.HLfit}}. 

The verbosity of the fitting functions is controlled by the \code{verbose["TRACE"]} value. 
The default is \code{TRACE=FALSE} (or 0) which is self-explanatory. \code{TRACE=TRUE} (or 1) shows values of outer-estimated parameters (and possibly fixed values of parameters that would be outer-estimated), some cryptic progress bar, and the attained value of the likelihood objective function (see \code{\link{inner-vs-outer}} for these two types of estimates). See Details for other values.

Tracing the progress of residual-dispersion models is distinctly controlled by \code{verbose["phifit"]}: see \code{\link{phi-resid.model}}.

If the fitted model includes a residual-dispersion model, some tracing output for the latter may be confusingly intermingled with tracing output of the mean-response model. The Details below are valid only for the mean-response model.

\code{confint{., verbose=.}} will not affect the verbosity of the fits performed by \code{confint} (but will rather affect messages specific to the progress of the \code{confint} procedure), \bold{unless} \code{verbose} is a named vector. In particular, \code{verbose["TRACE"]} may affect the fit (but this has been implemented only for fit objects inheriting from class \code{\link{pois4mlogit}}). \code{verbose["confint"]} may then be used to control \code{confint}. Thus, for \code{pois4mlogit} fits, \cr
\code{verbose=FALSE} and \code{verbose=c(confint=FALSE)} will affect only \code{confint} per se, in an equivalent way\cr
\code{verbose=c(TRACE=FALSE)} will affect only the fits performed to find the intervals, and\cr
\code{verbose=c(TRACE=FALSE, confint=FALSE)} will affect both.\cr
Keep in mind that the verbosity of \code{pois4mlogit} fits is also controlled by a \code{progress} argument. 


}

\details{
0<\code{TRACE}<1 only shows the cryptic progress bar.\cr \code{TRACE=1} was explained above. \cr \code{TRACE=2} will further show information about the progress of Levenberg-Marquardt-like steps for linear-predictor coefficients.\cr \code{TRACE=3} will further show information about the progress of distinct Levenberg-Marquardt-like steps random effects given fixed-effect coefficients.\cr \code{TRACE=4} will further show cryptic information about which matrix computations are requested.\cr \code{TRACE=5} will further report (through a call to the base \code{trace} function) the \code{str(.)} of the results of such matrix computations.\cr \code{TRACE=6} will further pause between iterations of the reweighted least-squares algorithm, allowing a browser session to be called.

(assuming no changes since this was first defined)
}
