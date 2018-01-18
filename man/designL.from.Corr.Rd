\name{designL.from.Corr}
\alias{designL.from.Corr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Computation of \dQuote{square root} of correlation matrix
}
\description{
  This function is not usually directly called by users, but arguments may be passed to it through higher-level calls (see Examples).  
  For given correlation matrix  \bold{C}, it computes a \dQuote{design matrix} \bold{L} such that \bold{C} = \bold{L} * t(\bold{L}). 
  \code{t(chol(.))} (Cholesky factorization) is a fast method for this computation, but it is not robust numerically and may even return an error, 
  in which cases more robust methods (\code{eigen} or \code{svd}) are used. Matrix roots are not unique (for example, they are lower triangular for \code{t(chol(.))}, and symmetric for \code{svd(.)}. As matrix roots are used to simulate samples under the fitted model (in particular in the parametric bootstrap implemented in \code{fixedLRT}), this implies that for given seed of random numbers, these samples will differ with these different methods (although their distribution should be identical). 
}
\usage{
designL.from.Corr(m = NULL, symSVD = NULL, try.chol = TRUE, try.eigen = FALSE, 
                  threshold = 1e-06, SVDfix = 1/10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{
  The matrix which 'root' is to be computed. This argument is ignored if \code{symSVD} is provided.
}
  \item{symSVD}{
  A list representing the symmetric singular value decomposition of the matrix which 'root' is to be computed. Must have elements \code{$u}, a matrix of eigenvectors, and \code{$d}, a vector of eigenvalues. 
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
  \item{SVDfix}{ 
  A solution to failures of \code{svd}: see Details.
}
}
\value{
  The \dQuote{square root of the input matrix}. Its rows and columns are labelled according to the columns of the original matrix.
}
\details{
  The function may call \code{\link{svd}}, for singular value decomposition (SVD) of a matrix \bold{M}. \code{svd} may return \dQuote{error code 1 from Lapack routine 'dgesdd'} (cf. unhelpful discussions on R forums). This can be circumvented by computing the SVD of \eqn{(1-x)}\bold{I}\eqn{+x}\bold{M} and deducing the 
  singular values of \bold{M} in a trivial way. The \eqn{x} value to be used in this fix is provided by the \code{SVDfix} argument.
  
  \code{svd} errors have occurred for correlation matrices that were close to the identity matrix except for a few large non-diagonal elements. Such matrices may occur in particular for the Matérn correlation model with low \eqn{\nu}, high \eqn{\rho}, and if some samples are spatially close. Then, an alternative fix to the \code{svd} problem may be to restrict the \eqn{\nu} and/or \eqn{\rho} ranges, using the \code{lower} and \code{upper} arguments of \code{corrHLfit}, although one should make sure that this has no bearing on the inferences. 
  
}

\examples{
\dontrun{
## try.chol argument passed to designL.from.Corr 
## through the '...' argument of higher-level functions
## such as HLCor, corrHLfit, fixedLRT:
data("scotlip")
HLCor(cases~I(prop.ag/10) +adjacency(1|gridcode)+offset(log(expec)),
      ranPars=list(rho=0.174),adjMatrix=Nmatrix,family=poisson(),
      data=scotlip,try.chol=FALSE)
}

}
