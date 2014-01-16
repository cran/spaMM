\name{HLCor}
\alias{HLCor}
\title{
                    Fits a (spatially) correlated mixed model, for given correlation parameters 
}
\description{
A convenient interface for \code{\link{HLfit}}, constructing the correlation matrix of random effects from the arguments, 
then estimating fixed effects and dispersion parameters using \code{HLfit}.
}
\usage{
HLCor(formula,ranPars,data,distMatrix,uniqueGeo,adjMatrix,
      corrMatrix,verbose=rep(FALSE,3),...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
   A \code{predictor}, i.e. a formula with attributes (see \code{\link{Predictor}}), or possibly simply a simple \code{formula} if an offset is not required.
  }
  \item{ranPars}{
  A list of values for correlation parameters (some of which are mandatory), and possibly also dispersion parameters 
  (optional, but passed to HLfit if present). 
  By default, the \code{\link{Matern.corr}} model is assumed,
  where the correlation parameters are \code{rho} (scale parameter(s)), \code{nu} (smoothness parameter), and (optional) \code{Nugget}. 
  The \code{rho} parameter can itself be a vector with different values for different geographic coordinates.
  If \code{corr.model} is \code{"adjacency"}, the only correlation parameter is a scalar \code{rho}.
  The dispersion parameters are the random effect variance, \code{lambda}, and the residual variance, \code{phi}.
}
  \item{data}{
   The data frame to be analyzed.
  }
  \item{distMatrix}{
    A distance matrix between geographic locations, forwarded to \code{Matern.corr}
}
  \item{uniqueGeo}{
    A matrix of non-redundant geographic locations. This is useful if the \code{rho} parameter is a vector with different values for different geographic coordinates, in which case
    a scaled distance matrix has to be reconstructed from \code{uniqueGeo} and \code{rho}. 
}
  \item{adjMatrix}{
    An adjacency matrix, used if the spatial random effects are of the form\cr
    \code{y ~ adjacency(1|<location index>)}. 
    Its rows and columns must be ordered as increasing values of the levels of the geographic 
    location index specifying the spatial random effect. For example, if the model formula 
    is\cr 
    \code{y ~ adjacency(1|geo.loc)} and \code{<data>$geo.loc} is 2,4,3,1,... the first row/column of the matrix refers to the fourth row of the data. 
}
  \item{corrMatrix}{
    A arbitrary (valid) correlation matrix, which each row correspond to levels of a variable <stuff> for the random effect (1|<stuff>) given in the model formula.
    This allows to analyze non-spatial model by giving for example a matrix of genetic correlations. \code{corrMatrix} will determine the \code{Lmatrix} in a Predictor object, which may be modified by additional matrices taking into account e.g. block effects as described in \code{\link{Predictor}}.
}
  \item{verbose}{
    A vector of up to three booleans, the first of which controls printed output from HLCor (currently not very useful), 
    the other ones being used as the 'verbose' argument of the HLfit call.  
}
  \item{\dots}{
    Further parameters passed to \code{HLfit} or to \code{designL.from.corr}.
}
}
\details{
  The correlation matrix for random effects can be specified using either the formula terms (not functions!) \code{Matern(1|<...>)} or \code{adjacency(1|<...>)} using the spatial coordinates in \code{<...>} (see Examples), or a user-provided distance matrix \code{distMatrix} to which the Matérn function is applied, or a user-provided correlation matrix \code{corrMatrix}.
  
  \code{Matern} will construct a correlation matrix according to the Matérn correlation function (see \code{\link{Matern.corr}}). 
  \code{adjacency} will construct a correlation matrix according to a conditional autoregressive (CAR) model, using the \code{adjMatrix}. 
  In the latter case the correlation matrix has the form (\bold{I}\eqn{-\rho}\code{adjMatrix})\eqn{^{-1}}. 
   This was implemented for comparison purposes as it is widely used, but is not really recommended for irregular lattices
  (see Wall, 2004 and Martellosio, 2012 for some insights on the implications of autoregressive models). 
}
\references{
  Wall M.M. (2004) A close look at the spatial structure implied by
the CAR and SAR models: Journal of Statistical Planning and Inference 121: 311-324.

 Martellosio, F. (2012) The correlation structure of spatial autoregressions, Econometric Theory 28, 1373-1391.
}

\value{
  A list with elements
  \item{hlfit}{The return value of the \code{HLfit} call}  
  \item{call}{the HLCor call}  
  \item{HLCor.info}{A list with optional additional information. May include the names of the coordinates used to compute the spatial 
                    correlation matrix between unique geographic locations.}  
}

\seealso{
\code{\link{Matern.corr}}, \code{\link{HLfit}}, \code{\link{corrHLfit}}
% possibly the currently undocumented designL.from.Corr
}
\examples{

#### Examples with Matérn correlations
data(Loaloa)
HLCor(cbind(npos,ntot-npos)~elev1+elev2+elev3+elev4+maxNDVI1+seNDVI
              +Matern(1|longitude+latitude),data=Loaloa,
      family=binomial(),ranPars=list(nu=0.5,rho=1/0.7)) ## takes ~ 6s

\dontrun{
#### example with an adjacency matrix (autoregressive model)
data(scotlip) 
HLCor(cases~I(prop.ag/10) +adjacency(1|gridcode)+offset(log(scotlip$expec)),
      ranPars=list(rho=0.174),
      adjMatrix=Nmatrix,family=poisson(),data=scotlip)
}

}
\keyword{ model }
