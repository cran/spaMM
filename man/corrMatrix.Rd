\name{corrMatrix}
\alias{corrMatrix}
\alias{distMatrix}

\title{Using the \code{corrMatrix} and \code{distmatrix} arguments}

\description{
\code{corrMatrix} is a formal argument of \code{HLCor}, also handled by higher-level fitting functions such as \code{fitme}, which can be used if the model formula contains a term of the form \code{corrMatrix(1|}\var{...}\code{)}. It describes a correlation matrix, possibly as a \code{dist} half-matrix object. The \code{\link{covStruct}} argument can be used for the same purpose and is much more general, in particular allowing to specify several correlation matrices.

A co\emph{variance} matrix (whose diagonal elements are not necessarily 1) can also be passed through this argument, but then its format must be a full matrix, not a \code{dist} object (as \code{dist} objects are interpreted as correlation matrices by filling the diagonal with 1's). 

The way the rows and columns of the matrix are matched to the rows of the \code{data} depends on the nature of the grouping term \var{...} in \code{corrMatrix(1|}\var{...}\code{)}, and on the \code{dimnames} of the matrix, both of which should be carefully controlled (see Examples).

\code{corrMatrix} specify a fixed correlation matrix. It can be used to fit parametric correlation models not already implemented in \pkg{spaMM}, by wrapping the call of the fitting function with given correlation matrix in an objective function called by a general optimizer such as \code{optim}. However, for the same purpose, it may be more efficient to use
a \code{\link{corrFamily}} \dQuote{constructor}. 

For geospatial models, the alternative \code{distMatrix} argument can also be used: a spatial correlation matrix will be recomputed from it and from spatial correlation parameters, fixed or not. This allows in principle to use distance matrices not derived from geographical coordinates. But not all matrices obtained by, say, applying the Matern correlation function to an arbitrary distance matrix are valid correlation matrices.     
}
\details{
The simplest case is illustrated in the first examples below: the grouping term is identical to a single variable which is present in the \code{data}, whose levels match the rownames of the \code{corrMatrix}. As illustrated, the order of the data does not matter in that case, because the factor levels are used to match the \code{data} rows to the appropriate row and columns of the \code{corrMatrix}. The \code{corrMatrix} may even contain rows (and columns) in excess of the levels of the grouping term, in  which case these rows are ignored. 

These convenient properties no longer hold when the grouping term is not a single variable from the \code{data} (see \code{nofactor} fit in Examples), or when its levels do not correspond to row names of the matrix. In these cases, (1) no attempt is made to match the \code{data} rows to the row and column names of the \code{corrMatrix} (such attempt could succeed only if the user had given names to the matrix matching those that the called function could create from the information in the \code{data}, in which case the user should find easier to specify a single variable that can be matched); (2) the order of \code{data} and \code{corrMatrix} matter: internally, a single factor variable is constructed from all levels of the variables in the grouping term (i.e., from all levels of \code{latitude} and \code{longitude}, in the \code{nofactor} fit), with levels 1,2,3... that are matched to rows 1,2,3... of the \code{corrMatrix}. Thus the first row of the data is always associated to the first row of the matrix; (3) further, the dimension of the matrix must match the number of levels implied by the grouping term. For example, one might consider the case of 14 response values but of correlations between only 7 levels of a random effect, with two responses for each level. Then the matrix must be of dimension 7x7.
}
\examples{
# These examples demonstrate correct use of the corrMatrix argument 
# by using it to fit a Matern spatial model and comparing the results 
# to the result using the alternative basic Matern(1|.) syntax 

# The first series of examples uses a simplified syntax that works 
# only because each geographical location is represented only once 
# in the data. This condition is relaxed afterwards.

data("blackcap") 

# Fit to be reproduced by the different syntaxes: 
fitme(migStatus ~ means+ Matern(1|longitude+latitude),data=blackcap,
      fixed=list(nu=0.6285603,rho=0.0544659))

# Here we manually reconstruct the correlation matrix of this fit:
MLcorMat <- MaternCorr(proxy::dist(blackcap[,c("longitude","latitude")]),
                        nu=0.6285603,rho=0.0544659)
# (Crucially, rownames of 'blackcap' define the dimnames of 'MLcorMat')

# We create a factor ordered as geographical locations are ordered 
# in the data, hence matching (in this simplified example) 
# their order in the correlation matrix: 
blackcap$name <- as.factor(rownames(blackcap))                

fitme(migStatus ~ means+ corrMatrix(1|name),data=blackcap,
      corrMatrix=MLcorMat)

# Correct results after permutation of the matrix:
# (here the 'dist' object has to be converted to 'matrix' to allow permutation).
perm <- sample(14)
pmat <- proxy::as.matrix(MLcorMat, diag=1)[perm,perm] 
# Crucially, info about the permutation is provided 
# by the dimnames of 'pmat' which allow a correct match 
# to levels of the 'name' variable.

fitme(migStatus ~ means+ corrMatrix(1|name),data=blackcap,
      corrMatrix=pmat)

# It is possible, but more risky, not to use a factor whose levels 
# are specifically defined to match the dimnames of the matrix. 
# Note the message if we don't use such a factor:
(nofactor  <- fitme(migStatus ~ means+ corrMatrix(1|longitude+latitude),
                    data=blackcap, corrMatrix=MLcorMat))
      
#### Case with several samples in each location: ####
if (spaMM.getOption("example_maxtime")>0.7 &&
    requireNamespace("IsoriX", quietly = TRUE)) {
# [the fits are quite fast, but requireNamespace("IsoriX"...)
#  is too slow on CRAN servers].

data("GNIPDataDE", package = "IsoriX")

## Fit to be reproduced by the different syntaxes:
fitme(source_value ~ 1 + Matern(1|long+lat), data=GNIPDataDE, 
      fixed=list(nu=0.75, rho=0.2))

## Reconstruct the distance matrix of this fit:
dat <- GNIPDataDE
longlat  <- paste0(dat$long,":",dat$lat)
# Use unique() to get unique geographical coordinates *and* unique names:
distmat <- as.matrix(dist(unique(dat[,c("long","lat")])))
dimnames(distmat)  <- list(unique(longlat), unique(longlat))
# Define matching factor:
dat$ID <- factor(longlat)

## Various possible syntaxes for the fit (last one more general):  
# corrMatrix(1|.) term + general 'covStruct' argument (see help("covStruct"))
fitme(source_value ~ 1 + corrMatrix(1|ID), data=dat, 
      covStruct=list(corrMatrix=MaternCorr(distmat,rho = 0.2, nu=0.75)))

# corrMatrix(1|.) term + more ad hoc 'corrMatrix' argument
fitme(source_value ~ 1 + corrMatrix(1|ID), data=dat, 
      corrMatrix=MaternCorr(distmat,rho = 0.2, nu=0.75))
      
# Matern(1|.) term + even more ad hoc 'distMatrix' argument
fitme(source_value ~ 1 + Matern(1|ID), data=dat, 
      distMatrix=distmat,  fixed=list(rho = 0.2, nu=0.75))  
      
# 'distMatrix' in covStruct=list(distMatrix=distmat) would not work, 
# but a general syntax allowing multiple distance matrices is:

fitme(source_value ~ 1 + Matern(1|ID), data=dat, 
      distMatrix=list("1"=distmat), # more than one element might be specified
      # covStruct=list(...), can be used if required for other random effects
      fixed=list(corrPars=list("1"=list(rho = 0.2, nu=0.75)))
     )      
}

}