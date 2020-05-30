\name{covStruct}
\alias{covStruct}
\alias{Predictor}
\alias{as_precision}
%\alias{subset_prec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Specifying correlation structures}
\description{
\code{covStruct} is a formal argument of \code{HLCor}, also handled by \code{fitme} and \code{corrHLfit}, that allows one to specify the correlation structure for different types of random effects, It is an alternative to other ad hoc formal arguments such as \code{corrMatrix} or \code{adjMatrix}. It replaces the deprecated function \code{Predictor(...)} which has served as an interface for specifying the design matrices for random effects in early versions of \code{spaMM}. 

The main use of \code{covStruct} is to specify the correlation matrix of levels of a given random effect term, or its inverse (a precision matrix). Assuming that the design matrix of each random effect term follows the structure \bold{ZAL} described in \code{\link{random-effects}}, it is thus an indirect way of specifying the \dQuote{square root} \bold{L} of the correlation matrix. The optional \bold{A} factor can also be given by the optional \code{"AMatrices"} attribute of \code{covStruct}.      

\code{covStruct} is a \code{list} of matrices with names specifying the type of matrix considered: \code{covStruct=list(corrMatrix=<some matrix>)} or \code{covStruct=list(adjMatrix=<some matrix>)}, where the \dQuote{corrMatrix} or \dQuote{adjMatrix} labels are used to specify the type of information provided (accordingly, the names can be repeated: \code{covStruct=list(corrMatrix=<.>,corrMatrix=<.>)}). 
\cr
The covariance structure of a \code{corrMatrix(1|<grouping factor>)} formula term can be specified in two ways (see Examples): either by a correlation matrix factor (\code{covStruct=list(corrMatrix=<some matrix>)}), or by a precision matrix factor \bold{Q} such that the covariance factor is \eqn{\lambda}\bold{Q}\eqn{^{-1}}, using the type name \code{"precision"}: \code{covStruct=list(precision=<some matrix>)}. 
The function \code{as_precision} can be used to perform the conversion from correlation information to precision factor (using a crude solve() that may not always be efficient), but fitting functions may also perform such conversions automatically.

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
\code{structure(list(<some matrix>,types="corrMatrix"))}.
}
\seealso{\code{\link{pedigree}} for a type of applications where declaring a precision matrix is useful.}
\examples{
\dontrun{
data("blackcap") 
# a 'dist' object can be used to specify a corrMatrix:  
MLdistMat <- MaternCorr(proxy::dist(blackcap[,c("latitude","longitude")]),
                        nu=0.6285603,rho=0.0544659) # a 'dist' object!
blackcap$name <- as.factor(rownames(blackcap))     
fitme(migStatus ~ means + corrMatrix(1|name), data=blackcap,
      corrMatrix=MLdistMat)

#### Same result by different input and algorithm:
fitme(migStatus ~ means + corrMatrix(1|name), data=blackcap,
      covStruct=list(precision=as_precision(MLdistMat)))

# Manual version of the same:
as_mat <- proxy::as.matrix(MLdistMat, diag=1) 
prec_mat <- solve(as_mat) ## precision factor matrix
fitme(migStatus ~ means + corrMatrix(1|name), data=blackcap,
      covStruct=list(precision=prec_mat))

# Since no correlation parameter is estimated, 
# HLcor(., method="ML")  is here equivalent to fitme()
}
}