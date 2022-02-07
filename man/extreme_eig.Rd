\name{extreme_eig}
\alias{extreme_eig}
\alias{regularize}
\title{
Utilities for regularization of a matrix
}
\description{
\code{regularize} can be used to regularize (nearly-)singular correlation matrices. It may also be used to  regularize covariance matrices but will not keep their diagonal constant. Use on other types of matrices may give nonsense. The regularization corrects the diagonal of matrices with high condition number so that the condition number of a corrected matrix is the maximum value specified by \code{maxcondnum}. For that purpose, it needs the extreme eigenvalues of the matrix, by default provided by the function \code{extreme_eig}. Calls functions from \pkg{RSpectra} if available, and falls back on base functions otherwise.
}
\usage{
extreme_eig(M, symmetric, required = TRUE)
regularize(A, EEV=extreme_eig(A,symmetric=TRUE), maxcondnum=1e12)
}
\arguments{
  \item{M}{
Square matrix. Sparse matrices of class \code{d[s|g]CMatrix} (and some others too) are handled (some vagueness, as if it fails for some matrix types, an alternative function shoudl be easy to define based on this one as template.
}
  \item{A}{
Square matrix as \code{M}, assumed symmetric. 
}
  \item{symmetric}{
Whether the matrix is symmetric. Helpful to select efficient methods for this case if the matrix class does not implies its symmetry.
}
  \item{required}{
Whether the computation should be attempted independently of the size of the matrix.
}
  \item{EEV}{ Two extreme eigenvalue in the return format of \code{extreme_eig}
}
  \item{maxcondnum}{ Target condition number when regularization is performed}
}
\value{
 \code{extreme_eig} returns a vector of length 2, the largest and the smallest eigenvalues in this order.
 \code{regularize} returns a matrix, possibly in sparse format.
}
\examples{
H10 <- Matrix::Hilbert(10)
extreme_eig(H10,symmetric=TRUE) # ratio > 1e13
rH10 <- regularize(H10)
extreme_eig(rH10,symmetric=TRUE)  # ratio = 1e12
}
