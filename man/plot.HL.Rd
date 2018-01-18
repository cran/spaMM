\name{plot.HLfit}
\alias{plot}
\alias{plot.HLfit}
\title{
Model checking plots for mixed models
}
\description{
This function provides diagnostic plots for residual errors from the mean model and for random effects. Plots for the mean models are similar 
to those for GLMs, as described in Lee et al. 2006. Plots for residual errors consider 
the \emph{standardized} deviance residuals (Lee et al. 2006, p.52), and plots for random effects likewise consider standardized values, i.e. 
each random deviate divided by \eqn{\sqrt(1-q)} where \eqn{q} is the corresponding leverage for \eqn{\lambda}.
}
\usage{
\method{plot}{HLfit}(x, which = c("mean", "ranef"), 
      titles = list(
          meanmodel=list(outer="Mean model",devres="Deviance residuals", 
                         absdevres="|Deviance residuals|", resq="Residual quantiles", 
                         devreshist="Deviance residuals"), 
          ranef=list(outer="Random effects and leverages",qq="Random effects Q-Q plot", 
                     levphi=expression(paste("Leverages for ",phi)), 
          levlambda=expression(paste("Leverages for ",lambda))) 
        ), 
      control= list(), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  The return object of an HLCor / HLfit / corrHLfit call.
  }
  \item{which}{ A vector of keywords for different types of plots. 
    By default, two types of plots are presented on different devices: diagnostic plots for mean values, and diagnostic plots for random effects.
    Either one can be selected using this argument. Use keyword \code{"predict"} for a plot of predicted response against actual response.
  }
  \item{titles}{
    A list of the \code{main} (inner and outer) titles of the plots. See the default value for the format.
  }
  \item{control}{
    A list of default options for the plots. Defaults are \code{pch="+"} and \code{pcol="blue"} for points, and \code{lcol="red"} for curves.
  }
  \item{\dots}{
    Options passed from \code{plot.HLfit} to \code{par}.  
  }
}
\details{
  The standardized deviance residuals are defined as the deviance residuals divided by \eqn{\phi\sqrt(1-q)}, 
  where \eqn{q} is the corresponding leverage for \eqn{\phi}, and the deviance residuals are defined as for a GLM. The leverages are zero for ML methods. Otherwise, they depend on the fitting method used, as defined in the Details of \code{\link{HLfit}}. The \code{PQL} and \code{EQL-} method use leverages obtained as diagonal elements of the \dQuote{hat} matrix; more elaborate methods will introduce corrections for non-Gaussian response and for non-Gaussian random effects; and \dQuote{(.,1)} methods will add another correction taking into account the variation of the GLM weights in the logdet Hessian term of restricted likelihood. 
  
  In principle the deviance residuals for the mean model should have a nearly Gaussian distribution hence form a nearly straight line on a Q-Q plot. However this is (trivially) not so for well-specified (nearly-)binary response data nor even for well-specified Poisson response data with moderate expectations. Hence this plot is not so useful.  
}
\value{
  Returns the input object invisibly.
}
\references{
  Lee, Y., Nelder, J. A. and Pawitan, Y. (2006). Generalized linear models with random effects: unified analysis via
  h-likelihood. Chapman & Hall: London.
}
\examples{
## see example for data("scotlip")
}
\keyword{hplot}
