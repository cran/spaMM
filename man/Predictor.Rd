\name{Predictor}
\alias{Predictor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Interface for model formulas 
}
\description{
  \code{Predictor(...)} essentially returns the \code{list(...)} of its allowed arguments, 
  with some minimal syntax checking.
  It serves as a unified interface for the set of descriptors for a linear predictor. 
  In particular, the return value includes the design matrix for random effects. 
  Also, in e.g. \code{\link{glm}} the formula and the offset are distinct arguments, which 
  are here grouped as two elements of the predictor (though this is waiting for an improved input syntax).
}
\usage{
Predictor(formula, offset = NULL, LMatrix = NULL, ZLMatrix = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
  a \code{\link{formula}}, which can include fixed effects and random effects
}
  \item{offset}{
  an offset term, on the scale of the linear predictor. 
}
  \item{LMatrix}{
  The \dQuote{square root} of the correlation matrix beteen unique locations, see Details.
}
  \item{ZLMatrix}{
  The design matrix for random effects, see Details.
}
}
\details{
  In a spatial model a vector  of correlated random effects Lv can be constructed from uncorrelated ones, v,
  as \code{L \%*\% v} for some matrix L. L is a square matrix which dimension is the number of spatial locations.
  Several observations may be taken in the same location, and the Z matrix tells which element of 
  Lv affects each observation. dim(Z) is therefore (number of observations,number of locations). 
  The ZL product matrix is reconstructed from the L matrix and the model formula, unless it is given as argument.
}
\value{
  A list with the same elements (but potentially different values) as in the input.
}

\examples{
# to include an offset
data(scotlip)
pred.example <- Predictor(cases~I(prop.ag/10) +adjacency(1|gridcode),
                          offset=log(scotlip$expec))
HLCor(ranPars=list(rho=0.174),adjMatrix=Nmatrix,family=poisson(),
      data=scotlip,predictor=pred.example)
}
