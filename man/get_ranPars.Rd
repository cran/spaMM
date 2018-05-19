\name{get_ranPars}
\alias{get_ranPars}
\title{
More functions to extract various components of a fit
}
\description{
\code{get_ranPars} returns various subsets of random-effect parameters.
}
\usage{
get_ranPars(object, which=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A object of class \code{HLfit}, as returned by the fitting functions in \code{spaMM}.}
  \item{which}{For \code{get_ranPars}, the only non-null value is \code{"corrPars"}, to return correlation parameters of random effects.}
  \item{\dots}{Other arguments that may be needed by some method.}
}
\value{
\code{get_ranPars(.,which="corrPars")} returns a (possibly nested) list of correlation parameters (or NULL if there is no such parameter). Top-level elements correspond to the different random effects. The list has \code{"type"} attribute having the same nested-list structure and describing whether and how the parameters where fitted: \code{"fix"} means they where fixed, not fitted; \code{"var"} means they were fitted by \code{HLfit}'s specific algorithms; \code{"outer"} means they were fitted by a generic optimization method.
}
\examples{
data("wafers")
m1 <- HLfit(y ~X1+X2+(1|batch),
          resid.model = ~ 1 ,data=wafers,HLmethod="ML")
get_ranPars(m1,which="corrPars") # NULL since no correlated random effect          
}
