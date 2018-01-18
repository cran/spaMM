\name{LRT}
\alias{anova}
\alias{anova.HLfit}
\alias{LRT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Likelihood ratio test of fixed effects.
}
\description{
\code{LRT} performs a likelihood ratio (LR) test between two model fits, the ``full'' and the ``null'' model fits, 
currently differing only in their fixed effects. Parametric bootstrap p-values can be computed, either using the raw bootstrap distribution of the likelihood ratio, or a bootstrap estimate of the Bartlett correction of the LR statistic.
This function differ from \code{fixedLRT} in its arguments (model fits for \code{LRT}, but all arguments required to fit the models for \code{fixedLRT}), and in the format of its return value. The function will stop or return possibly incorrect results for models differing beyond their fixed effects. By conceptual drift, \code{anova} works as an alias for \code{LRT}.  
}
\usage{
\method{anova}{HLfit}(object, object2, ..., method="")
LRT(object, object2, boot.repl=0, nb_cores=NULL, ...)  
}
\arguments{
  \item{object,object2}{Two models fits being compared (their order does not matter).}
  \item{boot.repl}{the number of bootstrap replicates. }
  \item{nb_cores}{Number of cores to use for bootstrap computation. The default is \code{spaMM.getOption("nb_cores")}, and 1 if the latter is NULL. \code{nb_cores=1} which prevents the use of parallelisation procedures.}
  \item{method}{For development purposes, not documented. }
  \item{\dots}{
  Further arguments passed to or from other methods. 
}
}
\details{
  A raw bootstrap p-value can be computed from the simulated distribution as \code{(1+sum(t >= t0))/(N+1)} where \code{t0} is the original likelihood ratio, \code{t} the vector of bootstrap replicates and \code{N} its length. See Davison & Hinkley (1997, p. 141) for discussion of the adjustments in this formula. However, a computationally more economical use of the bootstrap is to provide a Bartlett correction for the likelihood ratio test in small samples. According to this correction, the mean value \eqn{m} of the likelihood ratio statistic under the null hypothesis is computed (here estimated by a parametric bootstrap) and the original LR statistic is multiplied by \eqn{n/m} where \eqn{n} is the number of degrees of freedom of the test. 
}
\value{
An object of class \code{fixedLRT}, actually a list with as-yet unstable format, but here with typical elements (depending on the options)
\item{fullfit}{the HLfit object for the full model;} 
\item{nullfit}{the HLfit object for the null model;}
\item{basicLRT}{A data frame including values of the likelihood ratio chi2 statistic, its degrees of freedom, and the p-value;}

and, if a bootstrap was performed: 
\item{rawBootLRT}{A data frame including values of the likelihood ratio chi2 statistic, its degrees of freedom, and the raw bootstrap p-value;}
\item{BartBootLRT}{A data frame including values of the Bartlett-corrected likelihood ratio chi2 statistic, its degrees of freedom, and its p-value;}
\item{bootInfo}{a list with the following elements:
 \describe{
\item{bootreps}{A table of fitted likelihoods for bootstrap replicates;  }
\item{meanbootLRT}{The mean likelihood ratio chi-square statistic for boostrap replicates;  }
 } %- describe
} %- item
} %- value
\references{
 Bartlett, M. S. (1937) Properties of sufficiency and statistical tests. Proceedings of the Royal Society (London) A 160: 268-282.

 Davison A.C., Hinkley D.V. (1997) Bootstrap methods and their applications. Cambridge Univ. Press, Cambridge, UK.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
See also \code{\link{fixedLRT}}.
}
\examples{
data("wafers")
## Gamma GLMM with log link
m1 <- HLfit(y ~X1+X2+X1*X3+X2*X3+I(X2^2)+(1|batch),family=Gamma(log),
          resid.model = ~ X3+I(X3^2) ,data=wafers,HLmethod="ML")
m2 <- update(m1,formula.= ~ . -I(X2^2))
anova(m1,m2)
}
\keyword{ htest }
