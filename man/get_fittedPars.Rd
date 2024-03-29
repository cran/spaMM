\name{get_ranPars}
\alias{get_ranPars}
\alias{get_fittedPars}
\alias{remove_from_parlist}
\title{
Operations on lists of parameters
}
\description{

\code{\link{get_fittedPars}} returns estimated parameters. 

\code{get_ranPars} returns various subsets of random-effect parameters (correlation or variance parameters), as controlled by its \code{which} argument. 
It is one of several extractors for fixed or estimated parameters of  different classes of parameters, for which a quick guide is

\code{get_ranPars}: for random-effect parameters, excluding residual dispersion (with a subtlety for \code{corrFamily} models: see Details);\cr
\code{\link{VarCorr}}: alternative extractor for random-effect (co)variance and optionally residual variance, in a data frame format;\cr
\code{\link{residVar}}: for residual variance parameters, family dispersion parameters, or information about residual variance models;\cr
\code{\link{get_residVar}}: alternative extractor of residual variances with different features inherited from \code{\link{get_predVar}};\cr
\code{\link{get_inits_from_fit}}: extracts estimated parameters from a fit, in a different format from \code{get_fittedPars}.\cr

\code{remove_from_parlist} removes elements from a list of parameters, and from its \code{type} attribute. 
}
\usage{
get_fittedPars(object, partial_rC="rm", phiPars=TRUE)
get_ranPars(object, which=NULL, verbose=TRUE, 
            lambda_names = "Group.Term", ...)
remove_from_parlist(parlist, removand=NULL, rm_names=names(unlist(removand)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of class \code{HLfit}, as returned by the fitting functions in \code{spaMM}.}
  \item{partial_rC}{Controls handling of partially-fixed random coefficients. The default as set by \code{"rm"} is to remove the fixed values as for other parameters. But alternative option \code{"keep"} will keep the fixed value, and \code{NA} will replace it by a NA.}
  \item{phiPars}{Boolean: whether to include the parameters of any residual-dispersion model for \eqn{phi} (se \code{\link{phi-resid.model}}) in the \code{rdisPars} element of the returned list.}
  \item{which}{NULL or character string. Use \code{which="corrPars"} to get the correlation parameters. Use \code{which="lambda"} to get variances. see Details for the meaning of this for heteroscedastic models, and Value for other possible \code{which} values.}
  \item{\dots}{Other arguments that may be needed by some method.}
  \item{verbose}{Boolean: Whether to print some notes.}
  \item{parlist}{A list of parameters. see Details.}
  \item{removand}{Optional. A list of parameters to be removed from \code{parlist}.}
  \item{rm_names}{Names of parameters to be removed from \code{parlist}. Mandatory if \code{removand} is not given.}
  \item{lambda_names}{By default the names of the \code{lambda} vector are built from the Group (RHS of random effect term of the for (LHS|RHS)) and Term (variable from LHS). By setting a non-default value of \code{lambda_names} the names will be integer indices of the random-effect term in the model formula (currently, for \code{which="ranef_var"} or NULL.}
}
\details{
For heteroscedastic random effects, such as conditional autoregressive models, the variance parameter \dQuote{lambda} refers to a common scaling coefficient. For other random-effect models, \dQuote{lambda} typically refers to the single variance parameter.

\code{remove_from_parlist} is designed to manipulate structured lists of parameters, such as a list with elements \code{phi}, \code{lambda}, and \code{corrPars}, the latter being itself a list structured as the return value of \code{get_ranPars(.,which="corrPars")}. \code{parlist} may have an attribute \code{type}, also with elements \code{phi}, \code{lambda}, and \code{corrPars}... If given, \code{removand} must have the same structure (but typically not all the elements of \code{parlist}); otherwise, \code{rm_names} must have elements which match names of \code{unlist(names(parlist))}. 

If a corrFamily parameter is fixed through the formula term, as in \code{ARp(1|time, p=3, fixed=c(p2=0))}, the fixed parameter is not considered a model parameter and \code{get_ranPars} will not extract it from the object. However, the parameter will be extracted if it has been fixed through \code{fitme}'s \code{fixed} argument rather than through the formula term (see example in \code{\link{ARp}}).
}
\value{

\code{get_fittedPars} returns a list of model parameters, with possible elements: \code{beta} (fixed-effect coefficients); \code{lambda}, \code{phi}, \code{ranCoefs} and \code{corrPars} (same meaning as in \code{\link{fixed}} parameters); \code{hyper}, for \code{\link{multIMRF}} models;  
the residual-dispersion parameters \code{beta_prec}, \code{NB_shape} and \code{COMP_nu} when they are single scalars; 
and \code{rdisPars} for more complex residual-dispersion parameters. See the specific \code{\link{resid.model}} and \code{\link{phi-resid.model}} documentations for the \code{rdisPars} format, dependent on the nature of the residual-dispersion parameter being modelized. Use \code{residVar(., which="fam_parm")} to extract the vector of fitted values of the dispersion parameter. 

\code{get_ranPars(.,which="corrPars")} returns a (possibly nested) list of correlation parameters (or NULL if there is no such parameter). Top-level elements correspond to the different random effects. The list has a \code{"type"} attribute having the same nested-list structure and describing whether and how the parameters where fitted: \code{"fix"} means they where fixed, not fitted; \code{"var"} means they were fitted by \code{HLfit}'s specific algorithms; \code{"outer"} means they were fitted by a generic optimization method. 

\code{get_ranPars(.,which="lambda")} returns a vector of variance values, one per random effect, including both fixed, \dQuote{outer}- and \dQuote{inner}-optimized ones. The variances of random-coefficients terms with correlation parameters are not included.

\code{get_ranPars(.,which="outer_lambda")} returns only \dQuote{outer}-optimized variance parameters, ignoring those fitted by \code{HLfit}'s specific algorithms. 

\code{get_ranPars(.,which=NULL)} (the default) is not fully defined. It returns a list including the results of \code{which="lambda"} and \code{which="corrPars"}, but possibly other elements too. 

\code{get_ranPars(.,which="fitted")} is designed to provide fitted parameters with respect to which an information matrix is to be calculated (using \pkg{numDeriv}. It excludes fixed values, and has no \code{type} attribute. 

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
See \code{\link{get_fittedPars}}, \code{\link{VarCorr}}, \code{\link{residVar}}, \code{\link{get_residVar}}, or \code{\link{get_inits_from_fit}} as described in the quick guide above. 
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
