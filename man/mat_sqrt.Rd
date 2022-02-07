\name{mat_sqrt}
\alias{mat_sqrt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Computation of \dQuote{square root} of symmetric positive definite matrix
}
\description{

  \code{mat_sqrt} is not usually directly called by users, but arguments may be passed to it through higher-level calls (see Examples).  
  For given matrix \bold{C}, it computes a factor \bold{L} such that \bold{C} = \bold{L} * t(\bold{L}), handling issues with nearly-singular matrices.   The default behavior is to try Cholesky factorization, and use \code{\link{eigen}} if it fails. 
Matrix roots are not unique (for example, they are lower triangular for \code{t(chol(.))}, and symmetric for \code{svd(.)}. As matrix roots are used to simulate samples under the fitted model (in particular in the parametric bootstrap implemented in \code{fixedLRT}), this implies that for given seed of random numbers, these samples will differ with these different methods (although their distribution should be identical). 
}
\usage{
mat_sqrt(m = NULL, symSVD = NULL, try.chol = TRUE, condnum=1e12)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{
  The matrix whose 'root' is to be computed. This argument is ignored if \code{symSVD} is provided.
}
  \item{symSVD}{
  A list representing the symmetric singular value decomposition of the matrix which 'root' is to be computed. Must have elements \code{$u}, a matrix of eigenvectors, and \code{$d}, a vector of eigenvalues. 
}
  \item{try.chol}{
  If \code{try.chol=TRUE}, the Cholesky factorization will be tried. 
}
  \item{condnum}{ (large) numeric value.  In the case \code{chol()} was tried and failed, the matrix is regularized so that its (matrix 2-norm) condition number is reduced to \code{condnum} (in version 3.10.0 this correction has been implemented more exactly than in previous versions).
}
}
\value{
  For non-NULL \code{m}, its matrix root, with rows and columns labelled according to the columns of the original matrix.
If \code{eigen} was used, the symmetric singular value decomposition (a list with members \code{u} (matrix of eigenvectors) and \code{d} (vector of eigenvalues)) is given as attribute. 
}

\examples{
\dontrun{
## try.chol argument passed to mat_sqrt 
## through the '...' argument of higher-level functions
## such as HLCor, corrHLfit, fixedLRT:
data("scotlip")
HLCor(cases~I(prop.ag/10) +adjacency(1|gridcode)+offset(log(expec)),
      ranPars=list(rho=0.174),adjMatrix=Nmatrix,family=poisson(),
      data=scotlip,try.chol=FALSE)
}

}
