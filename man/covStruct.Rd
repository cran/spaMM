\name{covStruct}
\alias{covStruct}
\alias{Predictor}
\alias{as_precision}
%\alias{subset_prec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Specifying correlation structures}
\description{
\code{covStruct} is a formal argument of \code{HLCor}, also handled by higher level fitting functions such as \code{fitme} and \code{corrHLfit}, that allows one to specify the correlation structure for different types of random effects. It is an alternative to other ad hoc formal arguments such as \code{corrMatrix} or \code{adjMatrix}. It replaces the long-deprecated function \code{Predictor(...)} which has served as an interface for specifying the design matrices for random effects in early versions of \pkg{spaMM}. 

The main use of \code{covStruct} is to specify the correlation matrix of levels of a given random effect term. Assuming that the design matrix of each random effect term follows the structure \bold{ZAL} described in \code{\link{random-effects}}, it is thus an indirect way of specifying a \dQuote{square root} \bold{L} of the correlation matrix. The optional \bold{A} factor can also be given by the optional \code{"AMatrices"} attribute of \code{covStruct}.   

The correlation structure of a random effect can be specified directly as a correlation matrix, or by its inverse (a \emph{precision} matrix \bold{Q} such that the covariance matrix is \eqn{\lambda}\bold{Q}\eqn{^{-1}}), or by an \code{\link{adjacency}} matrix. \code{covStruct} is a \code{list} of matrices with names specifying the type of information provided:\cr 
\code{covStruct=list(corrMatrix=<some correlation (or covariance) matrix>)} or \code{covStruct=list(adjMatrix=<some adjacency matrix>)} or \code{covStruct=list(precision=<some precision matrix>)}. 

When several matrices are provided, these names can be repeated, and \code{NULL} list members may be necessary, e.g.\cr 
\code{covStruct=list(corrMatrix=<.>,"2"=NULL,corrMatrix=<.>)})\cr
when correlations matrices are required only for the first and third random effect. 

The function \code{as_precision} can be used to perform the conversion from correlation information to precision factor in a controlled way (but using a numerical inversion that may not always be efficient). Fitting functions may also perform such a conversion automatically. \bold{Do not run} \code{as_precision} on a matrix which is already the inverse of a correlation or covariance matrix, as the result will no longer be a precision matrix.


\code{"AMatrices"} is a list of matrices. The names of elements of the list does not matter, but the \emph{i}th \code{A} matrix, and its row names, should match the \emph{i}th \bold{Z} matrix, and its column names. This implies that \code{NULL} list members may be necessary, as for the \code{covStruct} list.
}
\usage{
as_precision(corrMatrix, condnum=1e12)  
}
\arguments{
\item{corrMatrix}{Correlation matrix, specified as \code{matrix} or as \code{dist} object.}
\item{condnum}{Numeric: when a standard Cholesky factorization fails, the matrix is regularized so that the regularized matrix has this condition number (in version 3.10.0 this correction has been implemented more exactly than in previous versions).}
}
\value{
\code{as_precision} returns a list with additional class \code{precision} and with single element a symmetric matrix of class \code{dsCMatrix}.
}
\details{
\code{covStruct} can also be specified as a list with an optional \code{"types"} attribute, e.g.\cr
\code{structure(list(<some matrix>,types="corrMatrix"))}.
}
\seealso{\code{\link{Gryphon}} and \code{\link{pedigree}} for a type of applications where declaring a precision matrix is useful.}
\examples{
\dontrun{
data("blackcap") 
# a 'dist' object can be used to specify a corrMatrix:  
MLdistMat <- MaternCorr(proxy::dist(blackcap[,c("latitude","longitude")]),
                        nu=0.6285603,rho=0.0544659) # a 'dist' object!
blackcap$name <- as.factor(rownames(blackcap))     
fitme(migStatus ~ means + corrMatrix(1|name), data=blackcap,
      corrMatrix=MLdistMat)

#### Same result using precision matrix,
(by_sp <- fitme(migStatus ~ means + corrMatrix(1|name), data=blackcap,
      covStruct=list(precision=as_precision(MLdistMat))))
# This typically induced the use of a distinct fitting algorithm:
how(by_sp) # "...using sparse-precision method..."

# Less controlled specification of precision matrix, 
# but essentially same result again:
as_mat <- proxy::as.matrix(MLdistMat, diag=1) 
prec_mat <- solve(as_mat) ## precision factor matrix
fitme(migStatus ~ means + corrMatrix(1|name), data=blackcap,
      covStruct=list(precision=prec_mat))
# !!! but do *not* use list(precision=as.precision(prec_mat)) !!!
 
# Since no correlation parameter is estimated, 
# HLcor(., method="ML")  is here equivalent to fitme()
}
}