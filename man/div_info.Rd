\name{div_info}
\alias{div_info}
\title{
Information about numerical problems
}
\description{
This function displays information about parameter values for which some numerical problems have occurred. It is prelimnary, has yet no good example of application, but some warnings might suggest its use.
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
  (logLik(numpb <- fitme(y ~ 1 + AR1(1|time), data = d, 
                         family = binomial(link=cauchit), method = "PQL")))
  div_info(numpb) # High lambda is often part of the problem,
                  # here only for negative AR parameter:      
  fitme(y ~ 1 + AR1(1|time), data = d, family = binomial(link=cauchit), 
        method = "PQL", lower=list(ARphi=0))
  # Less successful fix (or cover-up):
  fitme(y ~ 1 + AR1(1|time), data = d, family = binomial(link=cauchit), 
        method = "PQL", upper=list(lambda=20))
}
}
