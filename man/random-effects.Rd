\name{random-effects}
\alias{AMatrices}
\alias{random-effects}
\title{
Structure of random effects
}
\description{
The structure of random-effect models adjustable by \pkg{spaMM} can generally be described by the following steps. 

First, independent and  identically distributed (iid) random effects \bold{u} are drawn from one of the following distributions:  \bold{Gaussian} with zero mean, unit variance, and identity link; \bold{Beta}-distributed, where \eqn{u} \code{~} \eqn{B(1/(2\lambda),1/(2\lambda))} with mean=1/2, and var\eqn{=\lambda/[4(1+\lambda)]}; and with logit link \code{v=logit(u)};
 \bold{Gamma}-distributed random effects, where \eqn{u}\code{ ~ Gamma(shape=}1+1/\eqn{\lambda},\bold{scale=}1/\eqn{\lambda}): see \code{\link{Gamma}} for allowed links and further details; and \bold{Inverse-Gamma}-distributed random effects, where \eqn{u} \code{~} inverse-Gamma(\code{shape=}1+1/\eqn{\lambda},\bold{rate=}1/\eqn{\lambda}): see \code{\link{inverse.Gamma}} for allowed links and further details. 

Second, a transformation \bold{v}\eqn{=f}(\bold{u}) is applied (this defines \bold{v} whose elements are still iid). By default, \bold{v}=\bold{u} for gaussian random effects but not necessarily for other distributions of random effects (see \code{\link{Gamma}}).

For discussion of these alternative distributions, see Lee and Nelder 2001 or Lee et al. 2006, p. 178-.

Third, correlated random effects are obtained as \bold{Mv}, where the matrix \bold{M} can describe spatial correlation between observed locations, block effects (or repeated observations in given locations), and possibly also correlations involving unobserved locations (as is often the case for autoregressive models). In most cases \bold{M} is determined from the model formula, but it can also be controlled by \link{covStruct} argument. \bold{M} takes the form \bold{ZL} or \bold{ZAL}, where \bold{Z} is determined from the model formula, the optional \bold{A} factor is given by the optional \code{"AMatrices"} attribute of argument \code{\link{covStruct}} (see Examples), and \bold{L} can be determined from the model formula or from \code{covStruct}. In particular:

\verb{ * }\bold{Z} is typically, though not always, an incidence matrix: its elements \eqn{z_{ij}} are 1 if the \eqn{i}th observation is affected by the \eqn{j}th  element of \code{ALb}, and zero otherwise. The specification of \bold{Z} is controlled by the left- and right-hand side of the \dQuote{bar} formula terms \code{(<LHS>|<RHS>)}, as further described in \code{\link{spaMM-package}}.\cr

\verb{ * }For spatial random effects, \bold{L} is typically the Cholesky \dQuote{square root} of a correlation matrix determined by the random effect specification (e.g., \code{Matern(...)}), or given by the \code{covStruct} argument. This may be meaningful only for Gaussian random effects. Coefficients for each level of a random-coefficient model can also be represented as \bold{Lv} where \bold{L} is the \dQuote{square root} of a correlation matrix.\cr

\verb{ * }If there is one response value par location, \bold{L} for a spatial random effect is thus a square matrix whose dimension is the number of observations. Alternatively, several observations may be taken in the same location, and a matrix \bold{Z} (automatically constructed) tells which element of \bold{Lv} affects each observation. The linear predictor then contains a term of the form \bold{ZLv}, where \code{dim(Z)} is (number of observations,number of locations). 

\verb{ * }in \code{\link{IMRF}} random effects (IMRF for Interpolated Markov Random Fields), the realized random effects in response locations are defined as linear combinations \bold{ALv} of random effects \bold{Lv} in distinct locations. In that case the dimension of \bold{L} is the number of such distinct locations, an automatically constructed \bold{A} matrix maps them to the observed locations, and \bold{Z} again maps them to possibly repeated observations in observed locations.
}
\references{
Lee, Y., Nelder, J. A. (2001)  Hierarchical generalised linear models: A
synthesis of generalised linear models, random-effect models and structured
dispersions. Biometrika 88, 987-1006.

Lee, Y., Nelder, J. A. and Pawitan, Y. (2006). Generalized linear models with random effects: unified analysis via
h-likelihood. Chapman & Hall: London.
}
\examples{
set.seed(123)
fac <- factor(sample(letters[1:3],100,replace=TRUE))
toydata <- data.frame(fac=fac,
                      y=c(a=-1,b=2,c=0)[fac]+rnorm(100))
ranef(fitme(y~0+(1|fac), data=toydata)) # 3 values

## Merge levels "b" and "c" by (the odd way of) using an 'A' matrix:

# There is no non-trivial covariance structure, but the A matrix 
# can still be provided through the 'covStruct' argument
amatrix <- matrix(c(1,0,
                 0,1,
                 0,1), byrow=TRUE,ncol=2)
amatrix <-as(amatrix,"CsparseMatrix")
covstruct <- structure(list(), 
                       AMatrices=list("1"=amatrix))
                       
toyfit <- fitme(y~0+(1|fac), data=toydata, 
                covStruct=covstruct, method="REML")
ranef(toyfit) # 2 values
get_ZALMatrix(toyfit) # two-column incidence matrix for random effect 
}