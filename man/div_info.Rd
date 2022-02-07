\name{div_info}
\alias{div_info}
\title{
Information about numerical problems
}
\description{
This experimental function displays information about parameter values for which some numerical problems have occurred. Some warnings suggest its use.
}
\usage{
div_info(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of class \code{HLfit}, as returned by the fitting functions in \code{spaMM}.}
  \item{\dots}{Currently not used}
}
\value{Used mainly for the side effects (printed output) but returns invisibly either a single parameter vector (if a single numerical problem occurred) or a matrix of parameter ranges, or NULL if there is no problem to report.
}
\examples{
if (spaMM.getOption("example_maxtime")>25) {
  set.seed(1L)
  d <- data.frame(y = rbinom(100, 1, 0.5), time = 1:100)
  numpb1 <- fitme(y ~ 1 + AR1(1|time), data = d, init=list(lambda=NA), 
                         family = binomial(link=cauchit), method = "PQL")
  div_info(numpb1) # High lambda is often part of the problem,
                   # here only for negative AR parameter, as shown by the following 'fix':      
  #                 
  fitme(y ~ 1 + AR1(1|time), data = d, family = binomial(link=cauchit), 
        method = "PQL", lower=list(ARphi=0), init=list(lambda=NA)) # no div_info
  #
  # Less successfull attempt to fix the problem:
  numpb2 <- fitme(y ~ 1 + AR1(1|time), data = d, family = binomial(link=cauchit), 
        method = "PQL", upper=list(lambda=20), init=list(lambda=NA))  
  div_info(numpb2) # problem for negative AR parameter again
}
}
