\name{GHK_oneside}
\alias{GHK_oneside}
\concept{multivariate normal distribution}
\concept{GHK method}
\concept{integral}

\title{ Compute GHK estimate of multivariate normal one-sided integral}
\description{
  \code{GHK_oneside} computes the GHK estimate to the integral of a
  multivariate normal density over an intersection of half planes defined by a set
  of truncation points.  
}
\usage{
GHK_oneside(L, trunpt, above, nrep=1000)
}
\arguments{
  \item{L}{ lower triangular Cholesky root of Covariance matrix }
  \item{trunpt}{ vector of truncation points}
  \item{above}{ vector of indicators for truncation above(1) or below(0) }
  \item{nrep}{ number of draws to use in GHK }
}
\value{
A list with two elements: 
 \describe{
\item{logInt}{ The estimated logarithm of the integral  }
\item{seInt}{ The estimated standard error of \code{logInt}}
 } %- describe
}
\details{
  See Train (2009) or Genz and Bretz (2009) for an introduction to the GHK (Geweke–Hajivassiliou–Keane) method. 

  The C++ code underlying \code{GHK_oneside} builds on code from Peter Rossi, as implemented in version 2.2-5 of the \code{bayesm} package.
  It differs notably in the following ways:
  (i) it returns the log of the integral (which is more appropriate for values of the integral below the numeric limits for real numbers);
  (ii) it returns an estimates of the standard error of the log integral;
  (iii) it does not allow for computation of multiple integrals in one call.
  
}
\references{ 
  Genz, A. and Bretz, F. (2009), Computation of Multivariate Normal and t Probabilities. Lecture Notes in Statistics, Vol. 195. Springer-Verlag, Heidelberg. 

  Peter Rossi. (2012). bayesm: Bayesian Inference for Marketing/Micro-econometrics. R package version 2.2-5.
  
  Train, Kenneth (2009) Discrete choice methods with simulation. Cambridge UP, 2nd ed.
}
\examples{
Sigma=matrix(c(1,.5,.5,1),ncol=2)
L=t(chol(Sigma))
trunpt=c(1,1)
above=c(1,1)
GHK_oneside(L,trunpt,above,100)
}
\keyword{ distribution }
