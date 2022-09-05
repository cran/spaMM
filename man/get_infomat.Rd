\name{get_infomat}
\alias{get_infomat}
\title{
Information matrix of fitted model
}
\description{
This wraps \code{numDeriv::hessian} to provide the information matrix (negative Hessian matrix of log-likelihood) from the result of a \code{fitme} fit. This is mainly a stub for later developments, and the source code should be consulted to understand its limitations.
}
\usage{
get_infomat(object, verbose=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Result of \code{fitme} call.
}
  \item{verbose}{
Boolean: whether to print (as a list) the estimates of the parameters for which the Hessian will be computed, and a (sort of) progress bar.
}
  \item{\dots}{Arguments passed to \code{numDeriv::hessian}.}
}
\value{A matrix.}
\examples{
# This can be slow!

if (spaMM.getOption("example_maxtime")>10) {
  
  data("blackcap")
  foo <- fitme(migStatus ~ means+ Matern(1|longitude+latitude),data=blackcap) 
  (infomat <- get_infomat(foo))

  # Residual dispersion parameters are not yet handled:
  data("wafers")
  hlfit <- fitme(y ~ X1+X2+X1*X3+X2*X3+I(X2^2)+(1|batch), family=Gamma(log),
                 resid.model = ~ X3+I(X3^2) ,data=wafers)
  (infomat <- get_infomat(hlfit)) 

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
