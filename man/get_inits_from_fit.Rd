\name{get_inits_from_fit}
\alias{get_inits_from_fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Initiate a fit from another fit 
}
\description{
\code{get_inits_from_fit} is an extractor of some fitted values from a fit in a convenient format to initiate a next fit.  
}
\usage{
get_inits_from_fit(from, template = NULL, to_fn = NULL, inner_lambdas=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{from}{
Fit object (inheriting from class \code{"HLfit"}) from which fitted values are taken.
}
  \item{template}{
Another fit object. Usage with a \code{template} fit object is suitable for refitting this object using fitted values from the \code{from} object as starting values.
}
  \item{to_fn}{
  NULL or character: the name of the function to be used the next fit. If NULL, taken from \code{template} (if available), else from \code{from}. It is meaningful to provide a \code{to_fn} distinct from the function used to fit a \code{template}. 
} 
\item{inner_lambdas}{ Boolean; 
  Whether the output should include estimates of the dispersion parameters estimated by the iterative methods implemented in \code{HLfit}.
}
}
\value{
A list with elements 
\item{init, init.corrHLfit}{(depending on the fitting function) giving initial values for outer-optimization;} 
\item{init.HLfit}{giving initial values for the iterative algorithms in \code{HLfit}. It is itself a list with possible elements:
\describe{
\item{\code{fixef}}{for the coefficients of the linear predictor, adjusted to the format of the coefficients of the linear predictor of the \code{template} object, if available;}
\item{\code{ranCoefs}}{random-coefficients parameters (if \bold{not} outer-optimized).}
}
}
}
\seealso{
\code{\link{get_ranPars}} and \code{\link{VarCorr}}.
}
\examples{
\dontrun{
data("blackcap")
(corrhlfit <- corrHLfit(migStatus ~ means+ Matern(1|longitude+latitude),data=blackcap,
                        method="ML")) 
inits <- get_inits_from_fit(corrhlfit, to_fn = "fitme")                    
(fitfit <- fitme(migStatus ~ means+ Matern(1|longitude+latitude),data=blackcap, 
                  init=inits[["init"]])) 
inits <- get_inits_from_fit(corrhlfit, template = fitfit)                    
fitme(migStatus ~ means+ Matern(1|longitude+latitude),data=blackcap, 
      init=inits[["init"]])
# In these examples, inits$init.HLfit is useless 
# as it is ignored when LMMs are fitted by fitme().
}
}
