\name{member_extractors}
\alias{member_extractors}
\alias{get_ZALMatrix}
\title{
Member functions to extract various components of a fit
}
\description{
These extractors are defined as member functions of a fitted object. This allows their returned value to be stored in the fitted object after the first call to the extractor. Therefore, further calls do not need to recompute the result.    
\code{get_ZALMatrix} extracts the design matrix for the random effects \eqn{v}. More member extractors are available (\code{get_w_h_coeffs}, \code{get_beta_cov}, \code{get_beta_w_cov}, \code{get_invColdoldList}, \code{get_logdispObject}, and \code{get_info_crits}) but are undocumented.
}
\usage{
# object$get_ZALMatrix(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{The return object of a call to \code{HLfit} or similar fitting function.}
}
\value{
A matrix
}
\examples{
\dontrun{
 data(blackcap)
 thefit <- fitme(migStatus ~ means+ Matern(1|latitude+longitude),data=blackcap,
                    method="ML",fixed=list(nu=0.63, rho=0.054)) 
 thefit$get_ZALMatrix(thefit)
}
}
