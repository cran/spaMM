\name{verbose}
\alias{verbose}

\title{Tracking progress of fits}

\description{This (partially) documents the usage of the \code{verbose} argument of the fitting functions, and more specifically of \code{verbose["TRACE"]} values. 

Default is \code{TRACE=FALSE} (or 0) which is self-explanatory. \code{TRACE=TRUE} (or 1) shows values of outer-estimated parameters (and possibly fixed values of parameters that would be outer-estimated), some cryptic progress bar, and the attained value of the likelihood objective function (but when there inner-estimated dispersion parameters, the output is more difficult to describe concisely). Other values have effect may change in later versions without notice see Details).

If the fitted model includes a residual-dispersion mdoel, some tracing output for the latter may be confusingly intermingled with tracing output of the main-response model. The Details are valid only for the main-response model.
}

\details{
0<\code{TRACE}<1 only shows the cryptic progress bar.\cr  \code{TRACE=2} will further show information about the progress of Levenberg-Marquardt-like steps for linear-predictor coefficients.\cr \code{TRACE=3} will further show information about the progress of distinct Levenberg-Marquardt-like steps random effects given fixed-effect coefficients.\cr \code{TRACE=4} will further show cryptic information about which matrix computations are requested.\cr \code{TRACE=5} will further report (through a call to the base \code{trace} function) the \code{str(.)} of the results of such matrix computations.\cr \code{TRACE=6} will further pause between iterations of the reweighted least-squares algorithm, allowing a browser session to be called.
}
