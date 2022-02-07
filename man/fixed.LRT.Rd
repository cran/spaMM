\name{fixedLRT}
\alias{fixedLRT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Likelihood ratio test of fixed effects.
}
\description{
\code{fixedLRT} performs a likelihood ratio (LR) test between two models, the \dQuote{full} and the \dQuote{null} models, 
currently differing only in their fixed effects. Parametric bootstrap p-values can be computed, either using the raw bootstrap distribution of the likelihood ratio, or a bootstrap estimate of the Bartlett correction of the LR statistic.
This function differs from \code{LRT} in its arguments (model fits for \code{LRT}, versus all arguments required to fit the models for \code{fixedLRT}), and in the format of its return value. 
}
\usage{
fixedLRT(null.formula, formula, data, method, HLmethod = method, 
         REMLformula = NULL, boot.repl=0, control="DEPRECATED",
         control.boot="DEPRECATED", fittingFunction, seed=NULL,
         resp_testfn = NULL, ...)  
}
\arguments{
  \item{null.formula}{
Either a \code{formula} (as in \code{glm}) or a \code{predictor} (see \code{Predictor}) for the null model.
}
  \item{formula}{
Either a \code{formula} or a \code{predictor} for the full model.
}
  \item{data}{A data frame containing the variables in the model.}
  \item{method}{
   A method to fit the full and null models. 
   See \code{\link{method}} information about such methods. 
   The two most meaningful values of \code{method} in \code{fixedLRT} calls are: 
   \code{'ML'} for an LRT based on ML fits (generally recommended); and
   \code{'PQL/L'} for an LRT based on PQL/L fits (recommended for spatial binary data). 
   
   Also feasible, but more tricky, and not really recommended (see Rousset and Ferdy, 2014), is \code{'REML'}. 
   This will perform an LRT based on two REML fits of the data, *both* of which use the  
   same conditional (or \dQuote{restricted}) likelihood of residuals for estimating dispersion parameters \eqn{\lambda} and \eqn{\phi} (see \code{REMLformula} argument). 
   Further, REML will not be effective on a given dispersion parameter if a non-trivial init.corrHLfit value is provided for this parameter.  
}
 \item{HLmethod}{
   Kept for back-compatibility. Same as \code{method}, but may work only for \code{fittingFunction=corrHLfit}.
 }
  \item{REMLformula}{
  a formula specifying the fixed effects which design matrix is used in the REML correction 
  for the estimation of dispersion parameters, if these are estimated by REML. 
  This formula is by default that for the *full* model. 
}
  \item{boot.repl}{
  the number of bootstrap replicates.
}
  \item{control}{Deprecated.}
  \item{control.boot}{Deprecated.}
\item{fittingFunction}{Character string giving the function used to fit each model: either \code{"corrHLfit"} or \code{"fitme"}. Default is \code{"corrHLfit"} for small data sets (fewer than 300 observations), and \code{"fitme"} otherwise, but this may change in future versions.}
\item{seed}{Passed to \code{\link{simulate.HLfit}}}
\item{resp_testfn}{See argument \code{resp_testfn} of \code{\link{spaMM_boot}}}
%\item{type}{Passed to \code{\link{simulate.HLfit}}}
  \item{\dots}{
  Further arguments passed to or from other methods; presently, additional arguments passed to fitting functions.
}
}
\details{
  Comparison of REML fits is a priori not suitable for performing likelihood ratio tests. Nevertheless, it is possible to contrive them for testing purposes (Welham & Thompson 1997). This function generalizes some of Wehlam & Thompson's methods to GLMMs. 
  
  See Details in \code{\link{LRT}} for details of the bootstrap procedures.
}
\value{
An object of class \code{fixedLRT}, actually a list with as-yet unstable format, but here with  typical elements (depending on the options)
\item{fullfit}{the HLfit object for the full model;} 
\item{nullfit}{the HLfit object for the null model;}
\item{LRTori}{A likelihood ratio chi-square statistic} 
\item{LRTprof}{Another likelihood ratio chi-square statistic, after a profiling step, if any. }
\item{df}{the number of degrees of freedom of the test.}
\item{trace.info}{Information on various steps of the computation.  }

and, if a bootstrap was performed, the additional elements described in \code{\link{LRT}}. 
}

\references{
  Rousset F., Ferdy, J.-B. (2014) Testing environmental and genetic effects in the presence of spatial autocorrelation. Ecography, 37: 781-790.
    \doi{10.1111/ecog.00566}

  Welham, S. J., and Thompson, R. (1997) Likelihood ratio tests for fixed model
  terms using residual maximum likelihood, J. R. Stat. Soc. B 59, 701-714.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
See also \code{\link{corrHLfit}} and \code{\link{LRT}}.
}
\examples{
if (spaMM.getOption("example_maxtime")>1.9) {
 data("blackcap")
 ## result comparable to the corrHLfit examples based on blackcap
 fixedLRT(null.formula=migStatus ~ 1 + Matern(1|longitude+latitude),
       formula=migStatus ~ means + Matern(1|longitude+latitude), 
       method='ML',data=blackcap)
}
if (spaMM.getOption("example_maxtime")>156) {
 ## longer version with bootstrap
 fixedLRT(null.formula=migStatus ~ 1 + Matern(1|longitude+latitude),
       formula=migStatus ~ means + Matern(1|longitude+latitude), 
       method='ML',data=blackcap, boot.repl=100, seed=123) 
 }
}
\keyword{ htest }
