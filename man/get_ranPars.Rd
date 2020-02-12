\name{get_ranPars}
\alias{get_ranPars}
\alias{remove_from_parlist}
\title{
Operations on lists of parameters
}
\description{
\code{get_ranPars} returns various subsets of random-effect parameters.
\code{remove_from_parlist} removes elements from a list of parameters, and from its \code{type} attribute.
}
\usage{
get_ranPars(object, which=NULL, ...)
remove_from_parlist(parlist, removand=NULL, rm_names=names(unlist(removand)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of class \code{HLfit}, as returned by the fitting functions in \code{spaMM}.}
  \item{which}{For \code{get_ranPars}, the only non-null value is \code{"corrPars"}, to return correlation parameters of random effects.}
  \item{\dots}{Other arguments that may be needed by some method.}
  \item{parlist}{A list of parameters. see Details.}
  \item{removand}{Optional. A list of parameters to be removed from \code{parlist}.}
  \item{rm_names}{Names of parameters to be removed from \code{parlist}. Mandatory if \code{removand} is not given.}
}
\value{
\code{get_ranPars(.,which="corrPars")} returns a (possibly nested) list of correlation parameters (or NULL if there is no such parameter). Top-level elements correspond to the different random effects. The list has a \code{"type"} attribute having the same nested-list structure and describing whether and how the parameters where fitted: \code{"fix"} means they where fixed, not fitted; \code{"var"} means they were fitted by \code{HLfit}'s specific algorithms; \code{"outer"} means they were fitted by a generic optimization method. Note that \code{get_ranPars} does not provide random-coefficient parameters when these were fitted by \code{HLfit}'s specific algorithms, but that \code{\link{get_inits_from_fit}} can do that.

\code{remove_from_parlist} returns a list of model parameters with given elements removed, and likewise for its (optional) \code{type} attribute. See Details for context of application.
}
\details{
\code{remove_from_parlist} is designed to manipulate structured lists of parameters, such as a list with elements \code{phi}, \code{lambda}, and \code{corrPars}, the latter being itself a list structured as the return value of \code{get_ranPars(.,which="corrPars")}. \code{remove_from_parlist} may have an attribute \code{type}, also with elements \code{phi}, \code{lambda}, and \code{corrPars}... If given, \code{removand} must have the same structure (but typically not all the elements of \code{parlist}); otherwise, \code{rm_names} must have elements which match names of \code{unlist(names(parlist))}. 

The bias estimated by bootstrap is one of two terms included in the \code{get_cPredVar} computation, which together account for uncertainty in fixed effects and dispersion parameters. The computation does not account for uncertainty in correlation parameters \dQuote{outer-optimized} by \code{fitme} or \code{corrHLfit}. The bias computation would not be sufficient to account for the latter uncertainty, and thus the correlation parameters are fixed when the model is refitted on the bootstrap replicates. To account for uncertainty in correlation parameters, one should rather perform a parametric boostrap of the full model (typically using \code{spaMM_boot(., type="residual")}), which may take much more time.
}
\seealso{
\code{\link{get_inits_from_fit}}
}
\examples{
data("wafers")
m1 <- HLfit(y ~X1+X2+(1|batch),
          resid.model = ~ 1 ,data=wafers,HLmethod="ML")
get_ranPars(m1,which="corrPars") # NULL since no correlated random effect     

parlist1 <- list(lambda=1,phi=2,corrPars=list("1"=list(rho=3,nu=4),"2"=list(rho=5)))
parlist2 <- list(lambda=NA,corrPars=list("1"=list(rho=NA))) # values of elements do not matter
remove_from_parlist(parlist1,parlist2) ## same result as:
remove_from_parlist(parlist1,rm_names = names(unlist(parlist2)))
}
