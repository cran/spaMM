\name{lev2bool}
\alias{lev2bool}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Conversion of factor to 0/1 variable
}
\description{
It may be straightforward to add columns of indicator variables for each level of a factor to the data, by \cr 
\emph{<data>} \code{ <- cbind(}\emph{<data>}\code{, model.matrix( ~ }\emph{<factor>}\code{ - 1, data = }\emph{<data>}\code{))}.
Alternatively, indicator variables can be created on the fly for given levels, using the \code{lev2bool} function. 
}
\usage{
lev2bool(fac, lev)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fac}{
 An object coercible to \code{\link{factor}}.
}
  \item{lev}{
  The level of \code{fac} to be converted to 1.  
 }
}
\value{
A one-column matrix.
}
\examples{
## Elementary bivariate-response model

# Data preparation
#
fam <- rep(c(1,2),rep(6,2)) # define two biological 'families'
ID <- gl(6,2) # define 6 'individuals'
resp <- as.factor(rep(c("x","y"),6)) # distinguishes two responses per individual
set.seed(123)
toymv <- data.frame(
  fam = factor(fam), ID = ID, resp = resp, 
  y = 1 + (resp=="x") + rnorm(4)[2*(resp=="x")+fam] + rnorm(12)[6*(resp=="x")+as.integer(ID)]
)
toymv <- cbind(toymv, model.matrix( ~ resp - 1, data = toymv)) % also ade4::acm.disjonctif, maybe more memory-efficient.

# fit response-specific variances of random effect and residuals: 
#
(fitme(y ~ resp+ (0+respx|fam)+ (0+respy|fam), 
        resid.model = ~ 0+resp ,data=toymv))

# Same result by different syntaxes:

#  * by the lev2bool() specifier:
 (fitme(y ~ resp+ (0+lev2bool(resp,"x")|fam)+ (0+lev2bool(resp,"y")|fam), 
        resid.model = ~ 0+resp ,data=toymv))

#  * or by random-coefficient model using 'resp' factor:    
(fitme(y ~ resp+ (0+resp|fam), resid.model = ~ 0+resp ,data=toymv, 
       fixed=list(ranCoefs=list("1"=c(NA,0,NA)))))
       
#  * or by the dummy() specifier from lme4:
# (fitme(y ~ resp+ (0+dummy(resp,"x")|fam)+ (0+dummy(resp,"y")|fam), 
#        resid.model = ~ 0+resp ,data=toymv))
}
