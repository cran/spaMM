\name{get_ranPars}
\alias{get_ranPars}
\alias{remove_from_parlist}
\title{
Operations on lists of parameters
}
\description{
\code{get_ranPars} returns various subsets of random-effect parameters (correlation or variance parameters), as controlled by its \code{which} argument. A distinct documentation exists for \code{\link{get_inits_from_fit}} which is most useful to extract parameters from a fit in a form suitable to initialize another fit. \code{\link{VarCorr}} is yet another extractor for dispersion parameters.

\code{remove_from_parlist} removes elements from a list of parameters, and from its \code{type} attribute. 
}
\usage{
get_ranPars(object, which=NULL, ...)
remove_from_parlist(parlist, removand=NULL, rm_names=names(unlist(removand)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of class \code{HLfit}, as returned by the fitting functions in \code{spaMM}.}
  \item{which}{NULL or character string. Use \code{which="corrPars"} to get the correlation parameters. Use \code{which="lambda"} to get variances. see Details for the meaning of this for heteroscedastic models, and Value for other possible \code{which} values.}
  \item{\dots}{Other arguments that may be needed by some method.}
  \item{parlist}{A list of parameters. see Details.}
  \item{removand}{Optional. A list of parameters to be removed from \code{parlist}.}
  \item{rm_names}{Names of parameters to be removed from \code{parlist}. Mandatory if \code{removand} is not given.}
}
\details{
For heteroscedastic models, such as conditional autoregressive models, the variance parameter \dQuote{lambda} refers to a common scaling coefficient. For other random-effect models, \dQuote{lambda} typically refers to the single variance parameter.

\code{remove_from_parlist} is designed to manipulate structured lists of parameters, such as a list with elements \code{phi}, \code{lambda}, and \code{corrPars}, the latter being itself a list structured as the return value of \code{get_ranPars(.,which="corrPars")}. \code{parlist} may have an attribute \code{type}, also with elements \code{phi}, \code{lambda}, and \code{corrPars}... If given, \code{removand} must have the same structure (but typically not all the elements of \code{parlist}); otherwise, \code{rm_names} must have elements which match names of \code{unlist(names(parlist))}. 
}
\value{
\code{get_ranPars(.,which="corrPars")} returns a (possibly nested) list of correlation parameters (or NULL if there is no such parameter). Top-level elements correspond to the different random effects. The list has a \code{"type"} attribute having the same nested-list structure and describing whether and how the parameters where fitted: \code{"fix"} means they where fixed, not fitted; \code{"var"} means they were fitted by \code{HLfit}'s specific algorithms; \code{"outer"} means they were fitted by a generic optimization method. 

\code{get_ranPars(.,which="lambda")} returns a vector of variance values, one per random effect, including both \dQuote{outer}- and \dQuote{inner}-optimized ones.

\code{get_ranPars(.,which="outer_lambda")} returns only \dQuote{outer}-optimized variance parameters, ignoring those fitted by \code{HLfit}'s specific algorithms. 

\code{get_ranPars(. which="ranef_var")} (experimental) returns a list with elements
\describe{
    \item{\code{Var}}{same as \code{get_ranPars(.,which="lambda")}}
    \item{\code{lembda_est}}{A vector of variance values, one for each level of each random effect}
    \item{\code{outer}}{A vector or outer-optimized variance values, as returned by \code{get_ranPars(.,which="outer_lambda")}}
    \item{\ldots}{Other elements, subject to change in later versions.}
  }
\code{remove_from_parlist} returns a list of model parameters with given elements removed, and likewise for its (optional) \code{type} attribute. See Details for context of application.
}

\seealso{
\code{\link{get_inits_from_fit}} and \code{\link{VarCorr}}.
}
\examples{
data("wafers")
m1 <- HLfit(y ~X1+X2+(1|batch), resid.model = ~ 1, data=wafers, method="ML")
get_ranPars(m1,which="corrPars") # NULL since no correlated random effect     

parlist1 <- list(lambda=1,phi=2,corrPars=list("1"=list(rho=3,nu=4),"2"=list(rho=5)))
parlist2 <- list(lambda=NA,corrPars=list("1"=list(rho=NA))) # values of elements do not matter
remove_from_parlist(parlist1,parlist2) ## same result as:
remove_from_parlist(parlist1,rm_names = names(unlist(parlist2)))
}
