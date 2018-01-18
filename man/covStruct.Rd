\name{covStruct}
\alias{covStruct}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Specifying correlation structures}
\description{
\code{covStruct} is a formal argument of \code{HLCor}, also handled by \code{fitme} and \code{corrHLfit}, that allows one to specify the correlation structure for different types of random effects, It is an alternative to other ad hoc formal arguments such as \code{corrMatrix} or \code{adjMatrix}.\cr 
It is a \code{list} of matrices with names specifying the type of matrix considered: \code{covStruct=list(corrMatrix=<some matrix>)} or \code{covStruct=list(adjMatrix=<some matrix>)}, where the \dQuote{corrMatrix} or \dQuote{adjMatrix} names are used to specify the type of information provided.\cr
The covariance structure of a \code{corrMatrix(1|<grouping factor>)} formula term can be specified in two ways (see Examples): either by a correlation matrix factor (\code{covStruct=list(corrMatrix=<some matrix>)}), or by a precision matrix factor \bold{Q} such that the covariance is \eqn{\lambda}\bold{Q}\eqn{^{-1}}, using the type name \code{"precision"}: \code{covStruct=list(precision=<some matrix>)}. In this case, an algorithm efficient for \bold{sparse} precision matrices is used to fit the model.    
}
\details{
\code{covStruct} can also be specified as a list with an optional \code{"types"} attribute, e.g.\cr
\code{structure(list(<some matrix>,types="corrMatrix"))}\cr
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