\name{vcov}
\alias{beta_table}
\alias{vcov}
\alias{vcov.HLfit}
\alias{Corr}
\alias{VarCorr}
\alias{VarCorr.HLfit}
\title{
Extract covariance or correlation components from a fitted model object} 
\description{
\code{summary(<fit object>)$beta_table} returns the table of fixed-effect coefficients as it is printed by \code{summary}, including standard errors and t-values.

\code{vcov} returns the variance-covariance matrix of the fixed-effects coefficients (cf See Also for related computations involving random effects).

\code{Corr} by default returns correlation matrices of random effects (though see Details for user-defined correlation models).

\code{VarCorr} returns (co)variance parameters of random effects, and optionally the residual variance(s), from a fit object, in a data frame format roughly consistent with the method of objects of class \code{"lme"}, in particular including columns with consistent names for easier extraction. One may have to consult the summary of the object to check the meaning of the contents of this data frame (e.g., of 'variance' coefficients for non-gaussian random effects). Other extractors to consider are \code{\link{get_ranPars}} and \code{\link{get_inits_from_fit}}, the latter providing parameters in a form suitable for initializing a fit.

The covariance matrix of residuals of a fit can be obtained as a block of the hat matrix\cr 
(\code{\link{get_matrix}(., which="hat_matrix")}).  
This is (as other covariances matrices above) a matrix of expected values, generally assuming that the fitted model is correct and that its parameters are \dQuote{well} estimated, and should not to be confused with the computation of diagnostic correlations among inferred residuals of a fit.     
}
\usage{
\method{vcov}{HLfit}(object, ...)
\method{VarCorr}{HLfit}(x, sigma = 1, add_residVars=TRUE, verbose=TRUE, ...)
Corr(object, A=TRUE, cov2cor.=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object, x}{A fitted model object, inheriting from class \code{"HLfit"}, as returned by the fitting functions in \code{spaMM}.}
  %\item{which}{not documented, for development purposes}
  \item{add_residVars}{Boolean; whether to include residual variance information in the returned table.}
  \item{sigma}{ ignored argument, included for consistency with the generic function. }
  \item{A}{Boolean: Whether to return the correlation matrix described by the \bold{AL} matrix product, when there is an \bold{A} matrix (as for \code{IMRF} terms; see \code{\link{random-effects}}).}
  \item{cov2cor.}{Boolean: Whether to convert covariance matrices to correlation matrices (see Details).}
  \item{verbose}{Boolean: Whether to print some notes.}
  \item{\dots}{Other arguments that may be needed by some method.}
}
\details{
Any matrix returned by the \code{Corr} extractor is by default the unconditional correlation matrix of a vector \dQuote{\bold{ALv}} of random effects (as defined in  \code{\link{random-effects}}). 

But it may also be an heteroscedastic matrix for some random effects if \code{cov2cor.} is set to FALSE. In particular, the \code{\link{IMRF}} and \code{\link{MaternIMRFa}} models are by default defined in terms of the inverse of an heteroscedastic covariance matrix (with \code{tcrossprod} factor \bold{L}), and of a \bold{A} matrix of weights. The product \bold{AL} will be the \code{tcrossprod} factor of a covariance matrix rather than a correlation matrix, unless a non-default normalization was requested when declaring the random-effect terms. User-defined random-effects models may also be heteroscedastic. In all these cases \code{Corr} will by default return the correlation matrix, by applying \code{cov2cor} to the tcrossproduct. 
}
\value{
\code{vcov} returns a matrix. 

\code{Corr} returns a list, for the different random effect terms. For each random-effect term with nontrivial correlation structure, the returned element is a matrix, returned in base matrix format or in some class from \pkg{Matrix}. Otherwise the it is an information message.

\code{VarCorr} returns either NULL (if no variance to report, as for a poisson GLM) or a data frame with columns for the grouping factor, term, variance of random effect,  standard deviation (the root of the variance), and optionally for correlation of random effect in random-coefficient terms.  Information about the residual variance is optionally included as the last row(s) of the data frame, when relevant (gaussian- or Gamma-response models with single scalar parameter; beware the meaning of the residual variance parameter for Gamma-response models). 

Some variance parameters may be removed from the \code{VarCorr} output, with a message, such as the slope of the linear predictor describing the correlation model of an \code{adjacency} term (see \code{autoregressive}). The rare user of such parametrization should not consider this as a stable feature.   
}
\seealso{
\code{\link{get_inits_from_fit}} and \code{\link{get_ranPars}}. \cr

\code{\link{get_matrix}(., which="beta_v_cov")} for the joint covariance matrix of estimates/predictions of fixed-effect coefficients and random effects; and \cr
\code{\link{get_matrix}(., which="v_condcov")} for the covariance matrix of predictions of random effects given fixed effects (the variances corresponding to the \code{condsd} reported in some \pkg{lme4} output). Both of these computations refer to the random effects \bold{v} as defined in \code{\link{random-effects}}.
}
\examples{
data("wafers")
m1 <- HLfit(y ~ X1+X2+(1|batch), resid.model = ~ 1 ,data=wafers, method="ML")
vcov(m1)

# Example from VarCorr() documentation in 'nlme' package
data("Orthodont",package = "nlme")
sp1 <- fitme(distance ~ age+(age|Subject), data = Orthodont, method="REML")
VarCorr(sp1)
}
