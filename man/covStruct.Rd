\name{covStruct}
\alias{covStruct}
\alias{Predictor}
\alias{as_precision}
%\alias{subset_prec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Specifying correlation structures}
\description{
\code{covStruct} is a formal argument of \code{HLCor}, also handled by \code{fitme} and \code{corrHLfit}, that allows one to specify the correlation structure for different types of random effects, It is an alternative to other ad hoc formal arguments such as \code{corrMatrix} or \code{adjMatrix}. It replaces the deprecated function \code{Predictor(...)} which has served as an interface for specifying the design matrices for random effects in early versions of \code{spaMM}. 

It is assumed that the design matrices for the random effects take the form \bold{ZL} or \bold{ZAL} where the \bold{L} factor can be determined from the \code{covStruct} argument (or from the model formula; see Details), the \bold{Z} factor is determined from the model formula, and the optional \bold{A} factor is given by the optional \code{"AMatrices"} attribute of \code{covStruct}.

\code{covStruct} is a \code{list} of matrices with names specifying the type of matrix considered: \code{covStruct=list(corrMatrix=<some matrix>)} or \code{covStruct=list(adjMatrix=<some matrix>)}, where the \dQuote{corrMatrix} or \dQuote{adjMatrix} labels are used to specify the type of information provided (accordingly, the names can be repeated: \code{covStruct=list(corrMatrix=<.>,corrMatrix=<.>)}). 
\cr
The covariance structure of a \code{corrMatrix(1|<grouping factor>)} formula term can be specified in two ways (see Examples): either by a correlation matrix factor (\code{covStruct=list(corrMatrix=<some matrix>)}), or by a precision matrix factor \bold{Q} such that the covariance factor is \eqn{\lambda}\bold{Q}\eqn{^{-1}}, using the type name \code{"precision"}: \code{covStruct=list(precision=<some matrix>)}. In this case, an algorithm efficient for \bold{sparse} precision matrices is used to fit the model. 
The function \code{as_precision} can be used to perform the conversion from correlation information to precision factor (using a crude solve() that may not always be efficient).

\code{NULL} list members may be necessary, e.g.\cr 
\code{covStruct=list(corrMatrix=<.>,"2"=NULL,corrMatrix=<.>)})\cr
when correlations matrices are required only for the first and third random effect. 

\code{"AMatrices"} is a list of matrices. The names of elements of the list does not matter, but the \emph{i}th \code{A} matrix, and its row names, should match the \emph{i}th \bold{Z} matrix, and its column names. This implies that \code{NULL} list members may be necessary, as for the \code{covStruct} list.
}
\usage{
as_precision(corrMatrix)  
}
\arguments{
\item{corrMatrix}{Correlation matrix, specified as \code{matrix} or as \code{dist} object}
}
\value{
\code{as_precision} returns a list with addition class \code{precision} and with single element a symmetric matrix of class \code{dsCMatrix}.
}
\details{
\code{covStruct} can also be specified as a list with an optional \code{"types"} attribute, e.g.\cr
\code{structure(list(<some matrix>,types="corrMatrix"))}\cr


  In a spatial model, a vector of correlated random effects \bold{Lv} can be constructed from uncorrelated ones, \bold{v},
  for some matrix \bold{L} (this may be meaningful only for Gaussian random effects). Typically \bold{L} is the Cholesky \dQuote{square root} of a correlation matrix determined by the random effect specification (e.g., \code{Matern(...)}), or given by the \code{covStruct} argument. 
  
  If there is one response value par location, \bold{L} is thus a square matrix which dimension is the number of observations. Alternatively, several observations may be taken in the same location, and a matrix \bold{Z} automatically constructed) tells which element of \bold{Lv} affects each observation. The linear predictor then contains a term of the form \bold{ZLv}, where \code{dim(Z)} is (number of observations,number of locations). 
  
  Finally, in some applications the realized random effects in response locations may be viewed as linear combinations \bold{ALv} of random effects \bold{Lv} in distinct locations. In that case the dimension of \bold{L} is the number of such distinct locations, \bold{A} maps them to the observed locations, and \bold{Z} again maps them to possibly repeated observations in observed locations.
  
  Thus, in general the random term in the linear predictor is written \bold{Mv}, where  \bold{M=ZAL} is reconstructed from the element matrices (usually automatically constructed if needed).
  
 \code{attr(covStruct,"AMatrices")} should be used to specify \bold{A} matrices, and the other matrices are automatically constructed from the various arguments of a fit. 

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
HLCor(migStatus ~ means + corrMatrix(1|name), data=blackcap,
      covStruc=list(precision=as_precision(MLdistMat)), HLmethod="ML")
# Manual version of the same:
m <- proxy::as.matrix(MLdistMat, diag=1) 
prec_mat <- solve(m) ## precision factor matrix
HLCor(migStatus ~ means + corrMatrix(1|name), data=blackcap,
      covStruc=list(precision=prec_mat), HLmethod="ML")
}
}