\name{get_matrix}
\alias{get_matrix}
\alias{get_ZALMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Functions to extract various matrices used in fitting a model, from the resulting fit object. All augmented matrices follow (Henderson's) block order (upper blocks: X,Z; lower blocks: 0,I).
}
\description{
\code{get_matrix} is a first attempt at a unified extractor of various matrices from a fit.
\code{get_ZALMatrix} returns the design matrix for the random effects \eqn{v}.
}
\usage{
get_matrix(object, which="model.matrix", augmented=TRUE, ...)
get_ZALMatrix(object, as_matrix)
}
\arguments{
  \item{object}{A object of class \code{HLfit}, as returned by the fitting functions in \code{spaMM}.}
  \item{augmented}{Boolean; whether to return an augmented matrix for all model coefficients (fixed-effects coefficients and random-effect predictions) or only for fixed effects. Not operative for all \code{which} values (currently only for \code{which="left_ginv"}).}
  \item{which}{Which element to extract. For \code{"model.matrix"}, the design matrix for fixed effects (similarly to \code{stats::model.matrix}); for \code{"ZAL"}, the design matrix for random effects (same as code{get_ZALMatrix()}); for \code{"AugX"}, the (unweighted) augmented design matrix of the least-square problem; for \code{"hat_matrix"}, the projection matrix that gives model predictions from the (augmented) response vector; for \code{"left_ginv"}, the pseudo-inverse that gives the model coefficients from the (augmented) response vector. See Details for definitiosn and further options.}  
  \item{as_matrix}{Deprecated.}
  \item{\dots}{Other arguments that may be needed in some future versions of \code{spaMM}.}
}
\value{
A matrix, possibly in sparse-matrix format.
}
\details{
(Given the pain that is to to write maths in R documentation files, readers are gently asked to be tolerant about any imperfections of the following).

Model coefficents estimates of a (weighted) linear model can be written as \eqn{\bold{(X'WX)}^{-1}\bold{X'Wy}} where \bold{X} is the design matrix for fixed effects, \bold{W} a diagonal weight matrix, and \bold{y} the response vector. In a linear mixed model, the same expression holds in terms of Henderson's augmented design matrix, of augmented (still diagonal) weight matrix, and of augmented response vector. For GLMMs and hierarchical GLMs generally, the solution of each step of the iteratively reweighted least squares algorithm again has the same expression in terms of appropriately defined augmented matrices and vectors.  

\code{get_matrix} returns, for given values of the \code{which} argument, the following matrices from the model fit: \code{"AugX"}: \bold{X}; \code{"wei_AugX"}: \bold{WX}; \code{"wAugX"}: \eqn{\sqrt{\bold{W}}}\bold{X}; \code{"left_ginv"}: \eqn{\bold{X^{-}}=\bold{(X'WX)}^{-1}\bold{X'W}} (viewed as a pseudo-inverse since \\eqn{bold{X^{-}X}} is an identity matrix); \code{"hat_matrix"}: \eqn{\bold{X X^{-}}=\bold{X }\bold{(X'WX)}^{-1}\bold{X'W}}.


}