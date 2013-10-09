\name{fixedLRT}
\alias{fixedLRT}
%\alias{corrMM.LRT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Likelihood ratio test of fixed effects.
}
\description{
\code{fixedLRT} performs a likelihood ratio (LR) test between two models, the ``full'' and the ``null'' models, 
currently differing only in their fixed effects. This is the recommended interface to \code{corrHLfit} for performing such tests.
Both models are fit in same way, i.e. by ML, or allowing some REML correction for the estimation of dispersion parameters. 
A Bartlett correction of the LR statistic can be estimated using the bootstrap.

}
\usage{
fixedLRT(null.predictor,predictor,data,HLmethod,REMLformula=NULL,boot.repl=0,
          control=list(),control.boot=list(),...)  
}
\arguments{
  \item{null.predictor}{
Either a \code{formula} (as in \code{glm}) or a \code{predictor} (see \code{Predictor}) for the null model.
}
  \item{predictor}{
Either a \code{formula} or a \code{predictor} for the full model.
}
  \item{data}{A data frame containing the variables in the model.}
  \item{HLmethod}{
   A method to fit the full and null models. 
   See the identically-named \code{\link{HLfit}} argument for background information about such methods. 
   The two most meaningful values of \code{HLmethod} in \code{fixedLRT} calls are: 
   \code{'ML'} for an LRT based on ML fits (generally recommended); and
   \code{'PQL/L'} for an LRT based on PQL/L fits (recommended for spatial binary data). 

   
   Also feasible, but more tricky, and not really recommended (see Rousset and Ferdy, ms), is \code{'REML'}. 
   This will perform an LRT based on two REML fits of the data, *both* of which use the  
   same conditional (or \dQuote{restricted}) likelihood of residuals for estimating dispersion parameters \eqn{\lambda} and \eqn{\phi} (see \code{REMLformula} argument). 
   Further, REML will not be effective on a given dispersion parameter if a non-trivial init.corrHLfit value is provided for this parameter.  
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
  A set of control parameter for the fits of the data, mostly for development purposes. However, the ML fits may be slightly improved if initiated with the dispersion parameters estimates of an REML fit, in which case the option \code{control=list(REMLfits=TRUE)} may be useful (the effect appears small, however).   
}
  \item{control.boot}{
  Same as \code{control}, but for the fits of the bootstrap replicates. Again, the option \code{control.boot=list(REMLfits=TRUE)} may yield a small improvement in the fits, at the expense of more computation time.   
}
  \item{...}{
  Additional arguments passed to \code{corrHLfit}, including mandatory ones such as \code{data} and those ultimately passed to \code{designL.from.Corr}. With respect to the latter, note in particular that \code{try.chol} affects the simulation of samples for the parametric bootstrap, and although utlimate 
  differences in performance may be small, \code{try.chol=FALSE} may be slightly better.
}
}
\details{
  Comparison of REML fits is a priori not suitable for performing likelihood ratio tests. Nevertheless, it is possible to contrive them 
  for testing purposes (Wehlam & Thompson 1997). This function generalizes some of Wehlam & Thompson's methods to GLMMs. It further provides a Bartlett correction for the likelihood ratio test in small sample. According to this correction, the mean value \eqn{m} of the likelihood ratio statistic under the null hypothesis is computed (here estimated by a parametric bootstrap) and the original LR statistic is multiplied by \eqn{n/m} where \eqn{n} is the number of degrees of freedom of the test.
}
\value{
An object of class \code{fixedLRT}, actually a list with typical elements (depending on the options)
\item{fullfit}{the HLfit objects for the full model,} 
\item{nullfit}{the HLfit objects for the null model,} 
\item{LRTori}{A likelihood ratio chi-square statistic} 
\item{LRTprof}{Another likelihood ratio chi-square statistic, after a profiling step, if any. }
\item{df}{the number of degrees of freedom of the test.}
\item{trace.info}{Information on various steps of the computation.  }
\item{bootreps}{A table of fitted likelihoods for bootstrap replicates.  }
\item{meanbootLRT}{The mean likelihood ratio chi-square statistic for boostrap replicates.  }
}
\references{
  Bartlett, M. S. 1937. Properties of sufficiency and statistical tests. Proceedings of the Royal Society (London) A 160: 268-282.

  Rousset F., Ferdy, J.-B. (ms) Testing environmental and genetic effects in the presence of spatial autocorrelation.

  Welham, S. J., and Thompson, R. 1997. Likelihood ratio tests for fixed model
  terms using residual maximum likelihood, J. R. Stat. Soc. B 59, 701-714.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
See also \code{\link{corrHLfit}}.
}
\examples{
data(blackcap)
\dontrun{
## result comparable to the corrHLfit examples based on blackcap
fixedLRT(null.predictor=migStatus ~ 1 + Matern(1|latitude+longitude),
       predictor=migStatus ~ means + Matern(1|latitude+longitude), 
       HLmethod='ML',data=blackcap)

## longer version with bootstrap
fixedLRT(null.predictor=migStatus ~ 1 + Matern(1|latitude+longitude),
       predictor=migStatus ~ means + Matern(1|latitude+longitude), 
       HLmethod='ML',data=blackcap, boot.repl=100) 
}
}
\keyword{ htest }
