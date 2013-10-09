\name{plot.HLfit}
\alias{plot}
\alias{plot.HLfit}
\alias{plot.HLCor}
\alias{plot.corrHLfit}
\title{
Model checking plots for mixed models
}
\description{
This function provides diagnostic plots for residual errors from the mean model and for random effects. Plots for the mean models are similar 
to those for GLMs, as described in Lee et al. 2006. Plots for residual errors consider standardized residuals, that is, 
the standardized deviance residuals (Lee et al. 2006, p.52), and plots for random effects likewise consider standardized values, i.e. 
each random deviate divided by \eqn{\sqrt(1-q)} where \eqn{q} is the corresponding leverage for \eqn{\lambda}.
}
\usage{
\method{plot}{HLfit}(x, which = c(1, 2), 
      titles = list(
          meanmodel=list(outer="Mean model",devres="Deviance residuals", 
                         absdevres="|Deviance residuals|", resq="Residual quantiles", 
                         devreshist="Deviance residuals"), 
          ranef=list(outer="Random effects and leverages",qq="Random effects Q-Q plot", 
                     levphi=expression(paste("Leverages for ",phi)), 
          levlambda=expression(paste("Leverages for ",lambda))) 
        ), 
      control= list(), ...)
\method{plot}{HLCor}(x,...)
\method{plot}{corrHLfit}(x,...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  The return object of an HLCor / HLfit / corrHLfit call.
  }
  \item{which}{
    By default, two plots are presented. A subset can be selected using \code{which}.
  }
  \item{titles}{
    A list of the \code{main} (inner and outer) titles of the plots. See the default value for the format.
  }
  \item{control}{
    A list of default options for the plots. Defaults are \code{pch="+"} and \code{pcol="blue"} for points, and \code{lcol="red"} for curves.
  }
  \item{\dots}{
    Options passed from \code{plot.HLfit} to \code{par}, or from \code{plot.HLCor} or \code{plot.corrHLfit} to \code{plot.HLfit}.  
  }
}
\details{
  The standardized deviance residuals are defined as the deviance residuals divided by \eqn{\phi\sqrt(1-q)}, 
  where \eqn{q} is the corresponding leverage for \eqn{\phi}, and the deviance residuals are defined as for a GLM.
}
\note{
  This is currently the least thoroughly tested function of the package.
}
\value{
  Returns the input object invisibly.
}
\references{
  Lee, Y., Nelder, J. A. and Pawitan, Y. (2006). Generalised linear models with random effects: unified analysis via
  h-likelihood. Chapman & Hall: London.
}
\examples{
## see example for data(scotlip)
}
\keyword{hplot}