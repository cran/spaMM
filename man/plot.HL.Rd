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
      control= list(), ask=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class \code{HLfit}, as returned by the fitting functions in \code{spaMM}.}
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
  \item{ask}{
    Logical; passed to \code{devAskNewPage} which is run when a new device is opened by \code{code.HLfit}. 
  }
  \item{\dots}{
    Options passed from \code{plot.HLfit} to \code{par}.  
  }
}
\details{
  The standardized deviance residuals are defined as the deviance residuals divided by \eqn{\phi\sqrt(1-q)}, where the deviance residuals are defined as for a GLM, and the leverages considered here are those given by \code{hatvalues(., type="std")} (see \code{\link{hatvalues}} for details).
  
  In principle the deviance residuals for the mean model should have a nearly Gaussian distribution hence form a nearly straight line on a Q-Q plot. However this is (trivially) not so for well-specified (nearly-)binary response data nor even for well-specified Poisson response data with moderate expectations. Hence this plot is not so useful. The \code{DHARMa} package proposes better-behaved diagnostic plots (but the p-value that appears on one of these plots may not stand for a valid goodness-of-fit test). The current version of \code{DHARMa} should handle \code{spaMM} fit objects; otherwise, see \url{https://github.com/florianhartig/DHARMa/issues/95} for how to run \code{DHARMa} procedures on \code{spaMM} output.   
}
\value{
  Returns the input object invisibly.
}
\references{
  Lee, Y., Nelder, J. A. and Pawitan, Y. (2006). Generalized linear models with random effects: unified analysis via
  h-likelihood. Chapman & Hall: London.
}
\examples{
data("blackcap")
fit <- fitme(migStatus ~ 1+ Matern(1|longitude+latitude),data=blackcap,
             fixed=list(lambda=1,nu=1,rho=1))
plot(fit)
}
\keyword{hplot}
