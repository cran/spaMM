\name{sym_eigen}
\alias{sym_eigen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Singular Value Decomposition of a Symmetric Matrix
}
\description{
Computes the symmetric eigenvalue decomposition of a symmetric matrix, \eqn{X=U.D.U'}, where \eqn{U} is a matrix of orthogonal eigenvectors and \eqn{D} is a diagonal matrix of eigenvalues. The result is similar to that of \code{eigen(.,symmetric=TRUE)} but its evaluation uses the \code{SelfAdjointEigenSolver} function from the \code{Eigen} C++ library. No symmetry check is performed. 

This function was introduced to circumvent a bug in LAPACK (bug 113 in \url{http://www.netlib.org/lapack/bug_list.html}; originally \url{https://bugs.r-project.org/bugzilla3/show_bug.cgi?id=15211}). However, this bug has been fixed, and \code{sym_eigen} may be less accurate than the LAPACK one. It is no longer used is \code{spaMM} with default options (only the superseded function\code{\link{designL.from.Corr}}) may still call it). The structure of its return value was changed in spaMM version 2.4.123.  
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
  A list with members \code{vectors} (matrix of eigenvectors) and \code{values} (vector of eigenvalues).
}
\examples{
 hilbert <- function(n) { i <- 1:n; 1 / outer(i - 1, i, "+") }
 X <- hilbert(9)
 s <- sym_eigen(X)
 range(s$vectors \%*\% diag(s$values) \%*\% t(s$vectors) - X) # X=U D U'
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
