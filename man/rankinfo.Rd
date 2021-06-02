\name{rankinfo}
\alias{rankinfo}
\alias{get_rankinfo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Checking the rank of the fixed-effects design matrix
}
\description{
By default, fitting functions in \code{spaMM} check the rank of the design matrix for fixed effects, as \code{stats::lm}
or \code{stats::glm} do (but not, say, \code{nlme::lme}). This computation can be quite long.  
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
  \item{object}{An object of class \code{HLfit}, as returned by the fitting functions in \code{spaMM}.}
}
\details{
  The check is performed by a call to \code{qr()} methods for either dense or sparse matrices. 
  If the design matrix is singular, a set of columns from the design matrix that define a non-singular matrix is identified. Note that different sets may be identified by sparse- and dense-matrix \code{qr} methods.
}
\value{
  A list with elements \code{rank}, \code{whichcols} (a set of columns that define a non-singular matrix), and \code{method} (identifying the algorithm used).
}
\examples{
## Data preparation
# Singular matrix from ?Matrix::qr :
singX <- cbind(int = 1,
           b1=rep(1:0, each=3), b2=rep(0:1, each=3),
           c1=rep(c(1,0,0), 2), c2=rep(c(0,1,0), 2), c3=rep(c(0,0,1),2))
rownames(singX) <- paste0("r", seq_len(nrow(singX)))
donn <- as.data.frame(singX)
set.seed(123)
donn$y <- runif(6)

fitlm <- fitme(y~int+ b1+b2+c1+c2+c3,data=donn)
get_rankinfo(fitlm)

}