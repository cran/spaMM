\name{sym_eigen}
\alias{sym_eigen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Singular Value Decomposition of a Symmetric Matrix
}
\description{
Computes the symmetric eigenvalue decomposition of a symmetric matrix, \eqn{X=U.D.U'}, where \eqn{U} is a matrix of orthogonal eigenvectors and \eqn{D} is a diagonal matrix of eigenvalues. The result is similar to that of \code{eigen(.,symmetric=TRUE)} but its evaluation uses the \code{SelfAdjointEigenSolver} function from the \code{Eigen} C++ library. No symmetry check is performed. 
}
\usage{
sym_eigen(X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
  A symmetric matrix (possibly in \code{\link[Matrix]{sparseMatrix}} format).
}
}
\value{
  A list with members \code{u} (matrix of eigenvectors) and \code{d} (vector of eigenvalues).
}
\examples{
 hilbert <- function(n) { i <- 1:n; 1 / outer(i - 1, i, "+") }
 X <- hilbert(9)
 s <- sym_eigen(X)
 range(s$u \%*\% diag(s$d) \%*\% t(s$u) - X) # X=U D U'
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
