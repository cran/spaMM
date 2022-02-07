\name{random-effects}
\alias{random-effects}
\title{
Structure of random effects
}
\description{
The structure of random-effect models adjustable by spaMM can generally be described by the following steps. 

First, independent and  identically distributed (iid) random effects \bold{u} are drawn from one of the following distributions:  \bold{Gaussian} with zero mean, unit variance, and identity link; \bold{Beta}-distributed, where \eqn{u ~ B(1/(2\lambda),1/(2\lambda))} with mean=1/2, and var\eqn{=\lambda/[4(1+\lambda)]}; and with logit link \code{v=logit(u)};
 \bold{Gamma}-distributed random effects, where \eqn{u ~ }\code{Gamma(shape=}1+1/\eqn{\lambda},\bold{scale=}1/\eqn{\lambda}): see \code{\link{Gamma}} for allowed links and further details; and \bold{Inverse-Gamma}-distributed random effects, where \eqn{u ~ }inverse-Gamma(\code{shape=}1+1/\eqn{\lambda},\bold{rate=}1/\eqn{\lambda}): see \code{\link{inverse.Gamma}} for allowed links and further details.

Second, a transformation \bold{v}\eqn{=f}(\bold{u}) is applied (this defines \bold{v} whose elements are still iid). 

Third, correlated random effects are obtained as \bold{Mv}, where the matrix \bold{M} can describe spatial correlation between observed locations, block effects (or repeated observations in given locations), and correlations involving unobserved locations. In most cases \bold{M} is determined from the model formula, but it can also be controlled by \link{covStruct} argument. \bold{M} takes the form \bold{ZL} or \bold{ZAL}, where \bold{Z} is determined from the model formula, the optional \bold{A} factor is given by the optional \code{"AMatrices"} attribute of argument \code{\link{covStruct}} of \code{HLCor} (also handled by \code{fitme} and \code{corrHLfit}), and \bold{L} can be determined from the model formula or from \code{covStruct}. In particular:
\describe{
    \item{}{\bold{Z} is typically an incidence matrix: its elements \eqn{z_{ij}} are 1 if the \eqn{i}th observation is affected by the \eqn{j}th  element of \code{ALb}, and zero otherwise.}
    \item{}{For spatial random effects, \bold{L} is typically the Cholesky \dQuote{square root} of a correlation matrix determined by the random effect specification (e.g., \code{Matern(...)}), or given by the \code{covStruct} argument. This may be meaningful only for Gaussian random effects. Coefficients for each level of a random-coefficient model can also be represented as \bold{Lv} where \bold{L} is the \dQuote{square root} of a correlation matrix.}
    \item{}{If there is one response value par location, \bold{L} for a spatial random effect is thus a square matrix whose dimension is the number of observations. Alternatively, several observations may be taken in the same location, and a matrix \bold{Z} (automatically constructed) tells which element of \bold{Lv} affects each observation. The linear predictor then contains a term of the form \bold{ZLv}, where \code{dim(Z)} is (number of observations,number of locations). 
}
\item{}{in \code{\link{IMRF}} random effects (IMRF for Interpolated Markov Random Fields), the realized random effects in response locations are defined as linear combinations \bold{ALv} of random effects \bold{Lv} in distinct locations. In that case the dimension of \bold{L} is the number of such distinct locations, an automatically constructed \bold{A} matrix maps them to the observed locations, and \bold{Z} again maps them to possibly repeated observations in observed locations.
}
}
}
