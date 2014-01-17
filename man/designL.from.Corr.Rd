\name{designL.from.Corr}
\alias{designL.from.Corr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Computation of \dQuote{square root} of correlation matrix
}
\description{
  For given correlation matrix \bold{C}, this computes a \dQuote{design matrix} \bold{L} such that \bold{C} = \bold{L} * t(\bold{L}). 
  This function is normally not directly called by users, but they may need to pass options through higher-level functions (see Examples).
  \code{chol} (Cholesky factorization) is a fast method for this computation, but it is not robust numerically and may even return an error, 
  in which cases other methods (\code{eigen} or \code{svd}) are used. Note that these methods differ by much more than simply numerical accuracy,
  and in particular that this will affect the simulation of samples for the parametric bootstrap implemented in \code{fixedLRT}. 
}
\usage{
designL.from.Corr(m,try.chol=TRUE,try.eigen=FALSE,threshold=1e-06,debug=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{
  The matrix which 'root' is to be computed.
}
  \item{try.chol}{
  If \code{try.chol=TRUE}, the Cholesky factorization will be tried. 
}
  \item{try.eigen}{
  The default behavior is to try \code{\link{chol}}, and use \code{\link{svd}} if chol fails. 
  If \code{try.eigen=TRUE}, the \code{\link{eigen}} factorization will be tried before \code{svd}. 
  \code{eigen} is a compromise between speed and accuracy, but in our experience 
  it may *hang* so by default it is not tried.  
}
  \item{threshold}{ 
  A correction threshold for low eigenvalues is the case and eigensystem or singular-value decomposition are used.
}
  \item{debug}{ 
  Not documented, only for development purposes.
}
}
\value{
  The \dQuote{square root of the input matrix}. Its rows andcolumns are labelled according to the columns of the original matrix.
}


\examples{
## try.chol argument passed to designL.from.Corr 
## through the '...' argument of higher-level functions
## such as HLCor, corrHLfit, fixedLRT:

\dontrun{
#### using HLCor
data(scotlip)
HLCor(cases~I(prop.ag/10) +adjacency(1|gridcode)+offset(log(scotlip$expec)),
      ranPars=list(rho=0.174),adjMatrix=Nmatrix,family=poisson(),
      data=scotlip,try.chol=FALSE)
}

}
