\name{corrMatrix}
\alias{corrMatrix}

\title{Using a corrMatrix argument}

\description{
\code{corrMatrix} is an argument of \code{HLCor}, of calls \code{dist} or \code{matrix}, with is used if the model formula contains a term of the form \code{corrMatrix(1|<...>)}. It describes a correlation matrix, possibly as a \code{dist} object. A covariance matrix can actually be passed through this argument, but then it must be a full matrix, not a \code{dist} object. The way the rows and columns of the matrix are matched to the rows of the \code{data} depends on the nature of the grouping term {<...>}.
}
\details{
The simplest case is illustrated in the first two examples below: the grouping term is identical to a single variable which is present in the \code{data}, whose levels match the rownames of the \code{corrMatrix}. As illustrated by the second example, the order of the data does not matter in that case, because the factor levels are used to match the \code{data} rows to the appropriate row and columns of the \code{corrMatrix}. The \code{corrMatrix} may even contain rows (and columns) in excess of the levels of the grouping term, in  which case these rows are ignored. 

These convenient properties no longer hold when the grouping term is not a single variable from the \code{data} (third example below), or when its levels do not correspond to row names of the matrix. In these cases, (1) no attempt is made to match the \code{data} rows to the row and column names of the \code{corrMatrix}. Such attempt could succeed only if the user had given names to the matrix matching those that the the called function could create from the information in the \code{data}, in which case the user should find easier to specify a single variable that can be matched; (2) the order of \code{data} and \code{corrMatrix} matter; Internally, a single factor variable is constructed from all levels of the variables in the grouping term (i.e., from all levels of \code{latitude} and \code{longitude}, in the third example), with levels 1,2,3... that are matched to rows 1,2,3... of the \code{corrMatrix}. Thus the first row of the data is always associated to the first row of the matrix; (3) further, the dimension of the matrix must match the number of levels implied by the grouping term. For example, one might consider the case of 14 response values but of correlations between only 7 levels of a random effect, with two responses for each level. Then the matrix must be of dimension 7x7.
}
\examples{
data("blackcap") 
## Here we manually reconstruct the correlation matrix 
##  of the ML fit produced by corrHLfit:
MLcorMat <- MaternCorr(proxy::dist(blackcap[,c("latitude","longitude")]),
                        nu=0.6285603,rho=0.0544659)
blackcap$name <- as.factor(rownames(blackcap))                
## (1) Single variable present in the data 
HLCor(migStatus ~ means+ corrMatrix(1|name),data=blackcap,
      corrMatrix=MLcorMat,HLmethod="ML")
## (2) Same, permuted: still gives correct result
perm <- sample(14)
# Permuted matrix (with permuted names)
pmat <- as.matrix(MLcorMat)[perm,perm]
HLCor(migStatus ~ means+ corrMatrix(1|name),data=blackcap,
      corrMatrix=as.dist(pmat),HLmethod="ML")
## (3) Other grouping terms:
HLCor(migStatus ~ means+ corrMatrix(1|latitude+longitude),data=blackcap,
      corrMatrix=MLcorMat,HLmethod="ML")
}