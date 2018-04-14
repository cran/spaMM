\name{covStruct}
\alias{covStruct}
%\alias{subset_prec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Specifying correlation structures}
\description{
\code{covStruct} is a formal argument of \code{HLCor}, also handled by \code{fitme} and \code{corrHLfit}, that allows one to specify the correlation structure for different types of random effects, It is an alternative to other ad hoc formal arguments such as \code{corrMatrix} or \code{adjMatrix}.\cr 
It is a \code{list} of matrices with names specifying the type of matrix considered: \code{covStruct=list(corrMatrix=<some matrix>)} or \code{covStruct=list(adjMatrix=<some matrix>)}, where the \dQuote{corrMatrix} or \dQuote{adjMatrix} names are used to specify the type of information provided (accordingly, the names can be repeated: \code{covStruct=list(corrMatrix=<.>,corrMatrix=<.>)}).\cr
The covariance structure of a \code{corrMatrix(1|<grouping factor>)} formula term can be specified in two ways (see Examples): either by a correlation matrix factor (\code{covStruct=list(corrMatrix=<some matrix>)}), or by a precision matrix factor \bold{Q} such that the covariance factor is \eqn{\lambda}\bold{Q}\eqn{^{-1}}, using the type name \code{"precision"}: \code{covStruct=list(precision=<some matrix>)}. In this case, an algorithm efficient for \bold{sparse} precision matrices is used to fit the model.    

%Although the covariance matrix for a subset of levels is the corresponding subset of the covariance matrix, this is not true for precision matrices The \code{subset_prec} function can be used to perform the correct operations.
}
%\usage{
%subset_prec(prec_mat, corr_mat=NULL, rowcol_ids, tol = 1e-10, verbose = TRUE)  
%}
%\arguments{
%\item{prec_mat}{Precision matrix factor for the full set of levels of the random effect. Must be (Matrix::) sparse, and have dimnames.}
%\item{corr_mat}{Optional (see Details) correlation matrix factor for the full set of levels of the random effect. Must be in some Matrix:: format, and have dimnames.}
%\item{rowcol_ids}{Row and columns identifiers; must match the dimnames of the input matrix. Integer values are converted using \code{as.character()}}
%\item{tol}{A numerical tolerance value, passed to \code{Matrix::drop0} which is called after each matrix inversion operation.}
%\item{verbose}{Boolean: whether to print some information, or not.}
%}
\details{
\code{covStruct} can also be specified as a list with an optional \code{"types"} attribute, e.g.\cr
\code{structure(list(<some matrix>,types="corrMatrix"))}\cr
%\code{subset_prec} uses matrix inversions, which may be inefficient and exceed memory limits if the correlation matrix is numerically dense. If the \code{corr_mat} is not provided, it has to be computed by matrix inversion. It may be thus be useful to precompute the correlation matrix if \code{subset_prec} is called repeatedly on it.
}
\seealso{\code{\link{pedigree}} for a type of applications where declaring a precision matrix is useful.}
\examples{
\dontrun{
data("blackcap") 
# a 'dist' object can be used to specify a corrMatrix:  
MLdistMat <- MaternCorr(proxy::dist(blackcap[,c("latitude","longitude")]),
                        nu=0.6285603,rho=0.0544659) # a 'dist' object!
blackcap$name <- as.factor(rownames(blackcap))     
HLCor(migStatus ~ means + corrMatrix(1|name), data=blackcap,
      corrMatrix=MLdistMat, HLmethod="ML")
#### Same result by different input and algorithm:
m <- as.matrix(MLdistMat) ## leaves 0 on the diagonal! 
diag(m) <- 1 ## so that m is true correlation matrix 
prec_mat <- solve(m) ## precision factor matrix
HLCor(migStatus ~ means + corrMatrix(1|name), data=blackcap,
      covStruc=list(precision=prec_mat), HLmethod="ML")
}
}