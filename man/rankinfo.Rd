\name{rankinfo}
\alias{rankinfo}
\alias{get_rankinfo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Checking the rank of the fixed-effects design matrix
}
\description{
By default, fitting functions in \code{spaMM} check the rank of the design matrix for fixed effects, as \code{stats::lm}
or \code{stats::glm} do (but not, say, by \code{nlme::lme}). This computation can be quite long.  
To save time when fitting different models with the same fixed-effect terms to the same data, 
the result of the check can be extracted from a return object by \code{get_rankinfo()},
and can be provided as argument \code{control.HLfit$rankinfo} to another fit. Alternatively, the check will not be performed if
\code{control.HLfit$rankinfo} is set to \code{NA}.
}
\usage{
get_rankinfo(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A object of class \code{HLfit}, as returned by the fitting functions in \code{spaMM}.}
}
\details{
  The check is performed by a call to \code{qr()} methods for either dense or sparse matrices. 
  If the design matrix is singular, a set of columns from the design matrix that define a non-singular matrix is identified. Note that different sets may be identified by sparse- and dense-matrix \code{qr} methods.
}
\value{
  A list with elements \code{rank}, \code{whichcols} (a set of columns that define a non-singular matrix), and \code{method} (identifying the algorithm used).
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x) 
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
