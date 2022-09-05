\name{get_matrix}
\alias{get_matrix}
\alias{get_ZALMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract matrices from a fit
}
\description{
\code{get_matrix} is a first attempt at a unified extractor of various matrices from a fit. All augmented matrices follow (Henderson's) block order (upper blocks: X,Z; lower blocks: 0,I).
\code{get_ZALMatrix} returns the design matrix for the random effects \eqn{v}.
}
\usage{
get_matrix(object, which="model.matrix", augmented=TRUE, ...)
get_ZALMatrix(object, force_bind=TRUE)
}
\arguments{
  \item{object}{An object of class \code{HLfit}, as returned by the fitting functions in \code{spaMM}.}
  \item{augmented}{Boolean; whether to return a matrix for all model coefficients (augmented matrix for fixed-effects coefficients and random-effect predictions) or a matrix only for fixed effects. Not operative for all \code{which} values (currently only for \code{which="left_ginv"}).}
  \item{which}{Which element to extract. For \code{"model.matrix"}, the design matrix for fixed effects (similarly to \code{stats::model.matrix}); for \code{"ZAL"}, the design matrix for random effects (same as \code{get_ZALMatrix()}); for \code{"AugX"}, the (unweighted) augmented design matrix of the least-square problem; for \code{"hat_matrix"}, the projection matrix that gives model predictions from the (augmented) response vector; for \code{"left_ginv"}, the pseudo-inverse that gives the model coefficients from the (augmented) response vector. See Details for definitions and further options.}  
  \item{force_bind}{Boolean; with the non-default value \code{FALSE}, the function may return an object of class \code{\link{ZAXlist}}, which is poorly documented and for development purposes only. }
  \item{\dots}{Other arguments that may be needed in some future versions of \code{spaMM}.}
}
\value{
A matrix, possibly in \code{\link[Matrix]{sparseMatrix}} format.
}
\details{
(Given the pain that it is to write maths in R documentation files, readers are gently asked to be tolerant about any imperfections of the following).

Model coefficients estimates of a (weighted) linear model can be written as \bold{(X'WX)}\eqn{^{-1}}\bold{X'Wy} where \bold{X} is the design matrix for fixed effects, \bold{W} a diagonal weight matrix, and \bold{y} the response vector. In a linear mixed model, the same expression holds in terms of Henderson's augmented design matrix, of an augmented (still diagonal) weight matrix, and of an augmented response vector. For GLMMs and hierarchical GLMs generally, the solution of each step of the iteratively reweighted least squares algorithm again has the same expression in terms of appropriately defined augmented matrices and vectors.  

\code{get_matrix} returns, for given values of the \code{which} argument, the following matrices from the model fit:\cr 
\code{"AugX"}: \bold{X};\cr 
\code{"wei_AugX"}: \bold{WX};\cr 
\code{"wAugX"}: \eqn{\sqrt{}}{\sqrt}(\bold{W})\bold{X};\cr 
\code{"left_ginv"}: \bold{(X'WX)}\eqn{^{-1}}\bold{X'W} (the name stems from the fact that it is generalized inverse, denoted \bold{X}\eqn{^-}, since \bold{X}\bold{X}\eqn{^-}\bold{X}=\bold{X}, and it is a left one, since \bold{X}\eqn{^-}\bold{X} is an identity matrix when \bold{X} has full rank);\cr 
\code{"hat_matrix"}: \bold{XX}\eqn{^-}=\bold{X }\bold{(X'WX)}\eqn{^{-1}}\bold{X'W};\cr 
\code{"fixef_left_ginv"}: same as \code{"left_ginv"} but for the fixed-effect design matrix only (not to be confused with the corresponding block of \code{"left_ginv"}). 


}
