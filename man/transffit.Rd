\name{transffit}
\alias{transffit}
\alias{boxcox}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Transformation of the data.
}
\description{
\code{transffit} is an experimental helper function to fit a transformation of the data, 
such as the Box-Cox transformation which is the implemented default.
It is not particularly optimized, but is intended to help users to do it right.
When a non-default transformation is used, it is important to provide the corresponding \code{logDetJac} argument. 

}
\usage{
transffit(object, lower= -2, upper=2, verbose = FALSE, 
  extracts = quote({
    y <- object$y
  }), 
  updates = quote({
    if (lambda == 0) {
        newy <- log(y)
    } else newy <- (y^lambda - 1)/lambda
    refit <- update_resp(object, newresp = newy)
  }), 
  logDetJac = quote({
    (lambda - 1) * sum(log(y))
  }))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
 A fit object to the untransformed data, as produced by \code{fitme} or \code{fitmv}.
}
  \item{lower, upper}{
  Numeric vectors: optimization bounds (the function tentatively allows multiparameter transformations).
}
  \item{verbose}{
  Boolean: whether to print inforamtion about the progress of the procedure.
}
  \item{extracts}{
  A quoted R expression providing variables that should be available in the environment of the objective function, which refits the model on transformed variables. The default expression extracts the response variable so that it can be transformed by the objective function.
}
  \item{updates}{
  A quoted R expression providing the transformation (whose mandatory parameter names is \code{lambda}) applied to the response variable (and optionally to other variables), \emph{and} computes the refitted model (with mandatory named \code{refit}). 
}
  \item{logDetJac}{
  The logarithm of the determinant of the Jacobian of the transformation of the response variable, which is required to recover the likelihood, as the probability density of the original data, given the likelihod of the fit to the transformed data (Box and Cox 1964, eq.5). 
}
}
\value{
A list with elements the result of the optimization call, and the fit of the model to the best transformed data.
}
\references{
Box, G. E. P. and Cox, D. R. (1964) An analysis of transformations (with discussion). Journal of the Royal Statistical Society B, 26, 211–252.
}
\examples{
set.seed(123)
data(wafers)
wafers$x <- (rnorm(198, mean=(wafers$y)^(1/2),sd=0.1))^2
toyfit <- fitme(y ~x+(1|batch),data=wafers)

# Default transformation (Cox-Box)
transffit(toyfit, lower=-2,upper=2)

\dontrun{
# power transformation of response
transffit(toyfit, 
       updates = quote({   
         newy <- y^lambda
         refit <- update_resp(object, newresp = newy)
       }),
       logDetJac = quote({sum((lambda-1)*log(y)+log(lambda))}),
       lower=0.1, upper=2)

# Less standard power transformation of response *and* of predictor
# The update(object, data=.) syntax will be used
# so the 'data' must be extracted and updated:
transffit(toyfit,
       extracts = quote({
         locdata <- object$data
         y <- locdata$y
         x <- locdata$x
       }),
       updates = quote({   
         locdata$x <- x^lambda
         locdata$y <- y^lambda
         refit <- update(object, data=locdata)
       }),
       logDetJac = quote({sum((lambda-1)*log(y)+log(lambda))}),
       lower=0.1, upper=2)
}
}

