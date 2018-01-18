\name{fixedLRT}
\alias{fixedLRT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Likelihood ratio test of fixed effects.
}
\description{
\code{fixedLRT} performs a likelihood ratio (LR) test between two models, the ``full'' and the ``null'' models, 
currently differing only in their fixed effects. Parametric bootstrap p-values can be computed, either using the raw bootstrap distribution of the likelihood ratio, or a a bootstrap estimate of the Bartlett correction of the LR statistic.
This function differ from \code{LRT} in its arguments (model fits for \code{LRT}, but all arguments required to fit the models for \code{fixedLRT}), and in the format of its return value. 
}
\usage{
fixedLRT(null.formula, formula, data, 
         method, HLmethod = method, REMLformula=NULL, 
         boot.repl=0, control=list(),control.boot=list(), 
         fittingFunction, nb_cores=NULL, ...)  
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
   See \code{\link{HLfit}}'s \code{HLmethod} argument for background information about such methods. 
   The two most meaningful values of \code{method} in \code{fixedLRT} calls are: 
   \code{'ML'} for an LRT based on ML fits (generally recommended); and
   \code{'PQL/L'} for an LRT based on PQL/L fits (recommended for spatial binary data). 
   
   Also feasible, but more tricky, and not really recommended (see Rousset and Ferdy, 2014), is \code{'REML'}. 
   This will perform an LRT based on two REML fits of the data, *both* of which use the  
   same conditional (or \dQuote{restricted}) likelihood of residuals for estimating dispersion parameters \eqn{\lambda} and \eqn{\phi} (see \code{REMLformula} argument). 
   Further, REML will not be effective on a given dispersion parameter if a non-trivial init.corrHLfit value is provided for this parameter.  
}
 \item{HLmethod}{
   Kept for back-compatibility. Same as \code{method}, but wll work only for \code{fittingFunction=corrHLfit}.
 }
  \item{REMLformula}{
  a formula specifying the fixed effects which design matrix is used in the REML correction 
  for the estimation of dispersion parameters, if these are estimated by REML. 
  This formula is by default that for the *full* model. 
}
  \item{boot.repl}{
  the number of bootstrap replicates.
}
  \item{control}{
  A set of control parameters for the fits of the data, mostly for development purposes. However, if an initial value is provided for a dispersion parameter, a better one may be sought if further \code{control=list(prefits=TRUE)} (the effect appears small, however).   
}
  \item{control.boot}{
  Same as \code{control}, but for the fits of the bootstrap replicates. Again, the option \code{control.boot=list(prefits=TRUE)} may yield a small improvement in the fits, at the expense of more computation time.   
}
\item{fittingFunction}{Character string giving the function used to fit each model: either \code{"corrHLfit"} or \code{"fitme"}.}
\item{nb_cores}{Number of cores to use for bootstrap computation. The default is \code{spaMM.getOption("nb_cores")}, and 1 if the latter is NULL. \code{nb_cores=1} which prevents the use of parallelisation procedures.}
  \item{\dots}{
  Further arguments passed to or from other methods; in particular, additional arguments passed to \code{corrHLfit}, including mandatory ones such as \code{data} and those ultimately passed to \code{designL.from.Corr}. With respect to the latter, note that \code{try.chol} affects the simulation of samples for the parametric bootstrap, and although ultimate differences in performance may be small, \code{try.chol=FALSE} may be slightly better.
}
}
\details{
  Comparison of REML fits is a priori not suitable for performing likelihood ratio tests. Nevertheless, it is possible to contrive them for testing purposes (Wehlam & Thompson 1997). This function generalizes some of Wehlam & Thompson's methods to GLMMs. 
  
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
    \url{http://dx.doi.org/10.1111/ecog.00566}

  Welham, S. J., and Thompson, R. (1997) Likelihood ratio tests for fixed model
  terms using residual maximum likelihood, J. R. Stat. Soc. B 59, 701-714.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
See also \code{\link{corrHLfit}} and \code{\link{LRT}}.
}
\examples{
if (spaMM.getOption("example_maxtime")>2) {
 data("blackcap")
 ## result comparable to the corrHLfit examples based on blackcap
 fixedLRT(null.formula=migStatus ~ 1 + Matern(1|latitude+longitude),
       formula=migStatus ~ means + Matern(1|latitude+longitude), 
       HLmethod='ML',data=blackcap)
}
if (spaMM.getOption("example_maxtime")>223) {
 ## longer version with bootstrap
 fixedLRT(null.formula=migStatus ~ 1 + Matern(1|latitude+longitude),
       formula=migStatus ~ means + Matern(1|latitude+longitude), 
       HLmethod='ML',data=blackcap, boot.repl=100) 
 }
}
\keyword{ htest }
