\name{spaMM_glm.fit}
\alias{spaMM_glm.fit}
\alias{spaMM_glm}
\concept{regression}
\concept{logistic}
\concept{log-linear}
\concept{loglinear}
\title{Fitting generalized linear models without initial-value or divergence headaches}
\description{
\code{spaMM_glm.fit} is a stand-in replacement for \code{glm.fit}, which can be called through \code{glm} by using
\code{glm(<>, method="spaMM_glm.fit")}. Input and output structure are exactly as for \code{glm.fit}. It uses a Levenberg-Marquardt algorithm to prevent divergence of estimates. If the rcdd package is installed, the function can automatically find valid starting values or else indicate that no parameter value is feasible.
\code{spaMM_glm} is a convenient wrapper, calling \code{glm} with default method \code{glm.fit}, then calling method \code{spaMM_glm.fit}, with possibly different initial values, if \code{glm.fit} failed.
}
\usage{
spaMM_glm.fit(x, y, weights = rep(1, nobs), start = NULL, etastart = NULL, 
              mustart = NULL, offset = rep(0, nobs), family = gaussian(), 
              control = list(maxit=200), intercept = TRUE, singular.ok = TRUE)
spaMM_glm(formula, family = gaussian, data, weights, subset,
          na.action, start = NULL, etastart, mustart, offset,
          control = list(...), model = TRUE, method = c("glm.fit","spaMM_glm.fit"),
          x = FALSE, y = TRUE, singular.ok = TRUE, contrasts = NULL, strict=FALSE, ...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
All arguments except \code{strict} are common to these functions and their \code{stats} package equivalents, \code{glm} and \code{glm.fit}. Most arguments operate as for the latter functions, whose documentation is repeated  below. The \code{control} argument may operate differently. 

\item{formula}{an object of class \code{"\link{formula}"} (or one that
    can be coerced to that class): a symbolic description of the
    model to be fitted.  The details of model specification are given
    in the \sQuote{Details} section of \code{\link{glm}}.}
    
  \item{family}{a description of the error distribution and link
    function to be used in the model.  For \code{spaMM_glm} this can be a
    character string naming a family function, a family function or the
    result of a call to a family function.  For \code{spaMM_glm.fit} only the
    third option is supported.  (See \code{\link{family}} for details of
    family functions.)}
    
  \item{data}{an optional data frame, list or environment (or object
    coercible by \code{\link{as.data.frame}} to a data frame) containing
    the variables in the model.  If not found in \code{data}, the
    variables are taken from \code{environment(formula)},
    typically the environment from which \code{glm} is called.}  

  \item{weights}{an optional vector of \sQuote{prior weights} to be used
    in the fitting process.  Should be \code{NULL} or a numeric vector.}
    
  \item{subset}{an optional vector specifying a subset of observations
    to be used in the fitting process.}

  \item{na.action}{a function which indicates what should happen
    when the data contain \code{NA}s.  The default is set by
    the \code{na.action} setting of \code{\link{options}}, and is
    \code{\link{na.fail}} if that is unset.  The \sQuote{factory-fresh}
    default is \code{\link{na.omit}}.  Another possible value is
    \code{NULL}, no action.  Value \code{\link{na.exclude}} can be useful.}
    
  \item{start}{starting values for the parameters in the linear predictor.}

  \item{etastart}{starting values for the linear predictor.}

  \item{mustart}{starting values for the vector of means.}

  \item{offset}{this can be used to specify an \emph{a priori} known
    component to be included in the linear predictor during fitting.
    This should be \code{NULL} or a numeric vector of length equal to
    the number of cases.  One or more \code{\link{offset}} terms can be
    included in the formula instead or as well, and if more than one is
    specified their sum is used.  See \code{\link{model.offset}}.}

  \item{control}{a list of parameters for controlling the fitting
    process.  This is passed to \code{\link{glm.control}}, as for \code{glm.fit}. 
    Because one can assume that \code{spaMM_glm.fit} will converge in many cases where \code{glm.fit} does not, \code{spaMM_glm.fit} allows more iterations (200) by default. However, if \code{spaMM_glm.fit} is called through \code{glm(. . ., method="spaMM_glm.fit")}, then the number of iterations is controlled by the \code{glm.control} call within \code{glm}, so that it is 25 by default, overriding the \code{spaMM_glm.fit} default. } 
    
  \item{model}{a logical value indicating whether \emph{model frame}
    should be included as a component of the returned value.}
    
  \item{method}{A 2-elements vector specifying first the method to be used by \code{spaMM_glm} in the firt attempt at fitting the model, second the method to be used in a second attempt if the first failed. Possible methods include those shown in the default, \code{"model.frame"}, which returns the model frame and does no fitting, or user-supplied fitting functions. These functions can be supplied either as a function or a character string naming a function, with a function which takes the same arguments as \code{glm.fit}.  
  }
  

  \item{x, y}{For \code{spaMM_glm}: \code{x} is a design matrix of dimension
    \code{n * p}, and \code{y} is a vector of observations of length
    \code{n}.
    
    For \code{spaMM_glm.fit}: \code{x} is a design matrix of dimension
    \code{n * p}, and \code{y} is a vector of observations of length
    \code{n}.
  }
  
  \item{singular.ok}{logical; if \code{FALSE} a singular fit is an error.}
  
  \item{contrasts}{an optional list. See the \code{contrasts.arg}
    of \code{model.matrix.default}.}
    
  \item{intercept}{logical. Should an intercept be included in the
    \emph{null} model?}

  \item{strict}{logical. Whether to perform a fit by \code{spaMM_glm.fit} if \code{glm.fit} returned the warning \code{"glm.fit: algorithm did not converge"}.}
    
  \item{\dots}{
    arguments to be used to form the default
    \code{control} argument if it is not supplied directly.
  }  
  
}

\value{An object inheriting from class \code{glm}. See \code{\link[stats]{glm}} for details.}
\note{The source and documentation is derived in large part from those of \code{glm.fit}.}
\examples{
# glm() failures last checked in R-devel version (2017-04-19 r72556)
x <- c(8.752,20.27,24.71,32.88,27.27,19.09)
y <- c(5254,35.92,84.14,641.8,1.21,47.2)
\dontrun{ # fails:
  glm(y~ x,data=data.frame(x,y),family=Gamma(log)) 
}
spaMM_glm(y~ x,data=data.frame(x,y),family=Gamma(log))

## Gamma(inverse) examples
x <- c(43.6,46.5,21.7,18.6,17.3,16.7)
y <- c(2420,708,39.6,16.7,46.7,10.8)
\dontrun{ # fails (can't find starting value)
  glm(y~ x,data=data.frame(x,y),family=Gamma()) 
}
if (require("rcdd",quietly=TRUE)) {
  spaMM_glm(y~ x,data=data.frame(x,y),family=Gamma())
}
}
\keyword{models}
\keyword{regression}
