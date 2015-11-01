\name{HLCor}
\alias{HLCor}
\alias{corrMatrix}
\title{
                    Fits a (spatially) correlated mixed model, for given correlation parameters 
}
\description{
A convenient interface for \code{\link{HLfit}}, constructing the correlation matrix of random effects from the arguments, 
then estimating fixed effects and dispersion parameters using \code{HLfit}.
}
\usage{
HLCor(formula,ranPars=NULL,data,distMatrix,
      uniqueGeo=NULL,adjMatrix,
      corrMatrix,verbose=c(warn=TRUE,trace=FALSE,summary=FALSE),
      control.dist = list(),...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
   A \code{predictor}, i.e. a formula with attributes (see \code{\link{Predictor}}), or possibly simply a simple \code{formula} if an offset is not required.
  }
  \item{ranPars}{
  A list of values for correlation parameters (some of which are mandatory), and possibly also dispersion parameters 
  (optional, but passed to HLfit if present). 
  For the Matérn model, the correlation parameters are \code{rho} (scale parameter(s)), \code{nu} (smoothness parameter), and (optional) \code{Nugget}. 
  The \code{rho} parameter can itself be a vector with different values for different geographic coordinates.
  For the \code{"adjacency"} model, the only correlation parameter is a scalar \code{rho}.
  The dispersion parameters are the random effect variance, \code{lambda}, and the residual variance, \code{phi}.
}
  \item{data}{
   The data frame to be analyzed.
  }
  \item{distMatrix}{
    A distance matrix between geographic locations, forwarded to \code{MaternCorr}
}
  \item{uniqueGeo}{
    A matrix of non-redundant geographic locations. This is useful if the \code{rho} parameter is a vector with different values for different geographic coordinates, in which case
    a scaled distance matrix has to be reconstructed from \code{uniqueGeo} and \code{rho}. 
}
  \item{adjMatrix}{
    An adjacency matrix, used if a random effect of the form\cr
    \code{y ~ adjacency(1|<location index>)} is present. 
    Its rows and columns must be ordered as increasing values of the levels of the geographic 
    location index specifying the spatial random effect. For example, if the model formula 
    is\cr 
    \code{y ~ adjacency(1|geo.loc)} and \code{<data>$geo.loc} is 2,4,3,1,... the first row/column of the matrix refers to the fourth row of the data. 
}
  \item{corrMatrix}{
    A arbitrary (valid) correlation matrix, used if a random effect term of the form\cr
    \code{corrMatrix(1|<stuff>)} is present. Each row corresponds to levels of a variable <stuff>.
    This allows to analyze non-spatial model by giving for example a matrix of genetic correlations. Internally, \code{corrMatrix} will determine the \dQuote{\bold{L}} matrix, which may be modified by additional matrices (see Details in \code{\link{Predictor}}).
}
  \item{verbose}{
    A vector of booleans. \code{trace} controls various diagnostic (possibly messy) messages about the iterations.
    \code{summary} controls whether a summary of the fit is called by \code{HLfit}.
    \code{warn} is for programming purposes and best ignored.
}
  \item{control.dist}{
   A list of arguments that control the computation of the distance argument of the correlation functions. Possible elements are
   \describe{
   \item{rho.mapping}{ a set of indices controlling which elements of the \code{rho} scale vector scales which dimension(s) of the space in which (spatial) 
   correlation matrices of random effects are computed. See same argument in \code{\link{make_scaled_dist}} for details and examples.} 
   \item{dist.method}{\code{method} argument of \code{proxy::dist} function (by default, \code{"Euclidean"}, but e.g. \code{"Geodesic"} can be used to compute spherical distances.)}
   }
}
  \item{\dots}{
    Further parameters passed to \code{HLfit} or to \code{designL.from.corr}.
}
}
\details{
  The correlation matrix for random effects can be specified by various combination of formula terms and other arguments  (see Examples):
    \describe{
    \item{Basic Matérn model}{\code{Matern(1|<...>)} with \code{ranPars} argument, using the spatial coordinates in \code{<...>}. This will construct a correlation matrix according to the Matérn correlation function (see \code{\link{MaternCorr}});}
    \item{Matérn model with given distance matrix}{\code{Matern(1|<...>)} with \code{distMatrix} and \code{ranPars} argument;}
    \item{Given correlation matrix}{\code{corrMatrix(1|<...>)} with \code{corrMatrix} argument. In this case, a single factor variable is constructed from all levels of the joint variables in \code{<...>}. For example, one might consider the case of 14 response values but of correlations between only 7 levels of a random effect, with two responses for each level. Then the matrix is of dimension 7x7, and the indexing variable may be \code{idx=rep(seq(7),2)}, but also \code{dummy=rep(sample(7),2)}. If \code{dummy[1]=6}, then 6 is the first encountered level of the variable and it identifies the first row/column of the matrix, not the 6th. Then \code{dummy[8]=6} also means that the 8th response is affected by the same level of the random effect as the 1st response; }
    \item{CAR model with given adjacency matrix}{\code{adjacency(1|<...>)} with \code{adjMatrix} and \code{ranPars} argument. This will construct the 
    correlation matrix (\bold{I}\eqn{-\rho} \code{adjMatrix}\eqn{)^{-1}} for a conditional autoregressive (CAR) model. 
    This was implemented for comparison purposes as it is widely used, but is not really recommended for irregular lattices
    (see Wall, 2004 and Martellosio, 2012 for some insights on the implications of autoregressive models). See \code{\link{adjacency}} documentation for more information on fitting such models.}
  }
}
\references{
  Wall M.M. (2004) A close look at the spatial structure implied by
the CAR and SAR models: Journal of Statistical Planning and Inference 121: 311-324.

 Martellosio, F. (2012) The correlation structure of spatial autoregressions, Econometric Theory 28, 1373-1391.
}

\value{
  The return value of an \code{HLfit} call, with the following additional attributes: 
  \item{HLCorcall}{the HLCor call}  
  \item{info.uniqueGeo}{Unique geographic locations.}  
}

\seealso{
\code{\link{MaternCorr}}, \code{\link{HLfit}}, \code{\link{corrHLfit}}
% possibly the currently undocumented designL.from.Corr
}
\examples{

#### Matérn correlation using only the Matern() syntax
if (spaMM.getOption("example_maxtime")>2) {
 data(Loaloa)
 HLCor(cbind(npos,ntot-npos)~elev1+elev2+elev3+elev4+maxNDVI1+seNDVI
              +Matern(1|longitude+latitude),data=Loaloa,
      family=binomial(),ranPars=list(nu=0.5,rho=1/0.7)) 
}
%- : tested in simulate.HLCor.Rd

#### Using a corrMatrix (not necessarily Matérn)
data(blackcap) 
## Here we manually reconstruct the correlation matrix 
##  of the ML fit produced by corrHLfit:
MLcorMat <- MaternCorr(proxy::dist(blackcap[,c("latitude","longitude")]),
                        nu=0.6285603,rho=0.0544659)
# simple but easily misunderstood syntax:                        
HLCor(migStatus ~ means+ corrMatrix(1|latitude+longitude),data=blackcap,
      corrMatrix=MLcorMat,HLmethod="ML")
# An equivalent syntax worth understanding (see Details):      
bcdummy <- cbind(blackcap,dummy=sample(nrow(blackcap))) 
HLCor(migStatus ~ means+ corrMatrix(1|dummy),data=bcdummy,
      corrMatrix=MLcorMat,HLmethod="ML")
       
#### Matérn correlation using a distMatrix
MLdistMat <- as.matrix(proxy::dist(blackcap[,c("latitude","longitude")]))
HLCor(migStatus ~ means+ Matern(1|latitude+longitude),data=blackcap,
      distMatrix=MLdistMat,HLmethod="ML",ranPars=list(nu=0.6285603,rho=0.0544659))

#### example with an adjacency matrix (autoregressive model)
data(scotlip)
HLCor(cases~I(prop.ag/10) +adjacency(1|gridcode)+offset(log(scotlip$expec)),
      adjMatrix=Nmatrix,family=poisson(),data=scotlip)

}
\keyword{ model }
