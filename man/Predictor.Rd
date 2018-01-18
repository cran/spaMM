\name{Predictor}
\alias{Predictor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Interface for model formulas 
}
\description{
  \code{Predictor(...)} performs some minimal syntax checking, and returns a formula with attributes.
  It serves as a unified interface for the set of descriptors for a linear predictor, including design matrices for random effects. 
  
  In the current version this function could be ignored by users. 
}
\usage{
Predictor(formula, offset=NULL, LMatrix = NULL, AMatrix = NULL, ZALMatrix = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
  a \code{\link{formula}}, which can include fixed effects, random effects, and offsets.
}
  \item{offset}{
  a \code{\link{offset}} can be provided in this way, as a numeric vector. However, it \bold{may be better} to provide the offset as an 
  offset formula term (see \code{\link{scotlip}} example), in particular for later use with \code{\link{predict}} where the formula can be reevaluated on new data. 
}
  \item{LMatrix}{
  The \dQuote{square root} of the correlation matrix between unique locations, see Details.
}
  \item{AMatrix}{
  A matrix that relates observed (unique) locations to unobserved locations, see Details.
}
  \item{ZALMatrix}{Deprecated}
}
\details{
  In a spatial model a vector  of correlated random effects \bold{Lv} can be constructed from uncorrelated ones, \bold{v},
  for some matrix \bold{L} (this may be meaningful only for Gaussian random effects). Typically \bold{L} is the Cholesky \dQuote{square root} of a correlation matrix 
  determined by the random effect specification (e.g., \code{Matern(...)}), or given as the \code{corrMatrix} argument of \code{HLCor}. 
  
  If there is one realized random effect per response value, the linear predictor contains \bold{Lv}, where \bold{L} is a square matrix which dimension is the number of observations. 
  
  Several observations may be taken in the same location, and a matrix \bold{Z} (usually automatically constructed) tells which element of 
  \bold{Lv} affects each observation. The linear predictor then contains \bold{ZLv}, where \code{dim(Z)} is (number of observations,number of locations). 
  
  Finally, in some applications the realized random effects in response locations may be viewed as linear combinations \bold{ALv} of random effects \bold{Lv} in distinct locations. In that case the dimension of \bold{L} is the number of such distinct locations, \bold{A} maps them to the observed locations, and \bold{Z} again maps them to possibly repeated observations in observed locations.
  
  Thus, in general the random term in the linear predictor is written \bold{Mv}, where  \bold{M=ZAL} is reconstructed from the element matrices (usually automatically constructed if needed).
  
%- renvoyer a predictForNewX ...
}
\value{
  A formula with attributes. This return object has classes \code{formula} and \code{predictor}.
}

\examples{
# In the current version this function can be ignored by users,
# so examples are not required.
# (Use of AMatrix could perhaps be shown) 
}
