\name{plot.HLfit}
\alias{plot.HLfit}
\alias{`DHARMa-on-HLfit`}
\title{
Model checking plots for mixed models
}
\description{
This function provides diagnostic plots for residual errors from the mean model and for random effects. Plots for the mean models are similar to those for GLMs. However, \code{"std_dev_res"} ended as the type of the residuals plotted by the old version of plot.HLfit instead of the intended \code{"std_dev_rt"} (the standardized deviance residuals as described by Lee et al. 2006, p.52), by historical accident and lack of interest for the dubious results even with \code{"std_dev_rt"} (see Details). The \code{"RQR"} type may be more interesting type, for count data in particular.

Plots for random effects likewise consider standardized values.
}
\usage{
\method{plot}{HLfit}(x, 
    which=c("mean","ranef"),
    res_type="std_dev_rt",
    form = residuals(., type=res_type) ~ fitted(.),
    titles = list(
      meanmodel=list(outer="Mean model",devres="residuals", absdevres="|residuals|",
                     resq="Residual quantiles", devreshist="residuals"),
      ranef=list(outer="Random effects and leverages",qq="Random effects Q-Q plot",
                 levphi=expression(paste("Leverages for ",phi)), 
                 levlambda=expression(paste("Leverages for ",lambda)))
    ),
    control = list() , ask=TRUE, ...)
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
  \item{res_type}{
  Character string: the type of residual in the plots, passed to \code{\link{residuals.HLfit}(., type)}. 
  }
  \item{form}{
Either a formula specifying the desired type of plot for residuals, wherein 
%Any variable present in the original data frame used to obtain \code{x} can be referenced. In addition,
the fit object \code{x} can be referenced using the symbol \dQuote{.}; or (more for programming purposes) directly a vector of values for the residuals. 
%Conditional expressions on the right of a | operator can be used to define separate panels in a lattice display.
  }
  \item{\dots}{
    Options passed from \code{plot.HLfit} to \code{par}.  
  }
}
\details{
  In principle the standardized deviance residuals for the mean model should have a nearly Gaussian distribution hence form a nearly straight line on a Q-Q plot. However this is (trivially) not so for well-specified (nearly-)binary response data nor even for well-specified Poisson response data with moderate expectations. Hence this plot is not so useful. The \code{DHARMa} package proposes better-behaved diagnostic plots, but the p-value that appears on one of these plots may not stand for a valid goodness-of-fit test; see instead the \code{\link{gof}} procedure using the randomized quantile residuals (\code{"RQR"}). The current version of \code{DHARMa} should handle \code{spaMM} fit objects; otherwise, see \url{https://github.com/florianhartig/DHARMa/issues/95} for how to run \code{DHARMa} procedures on \code{spaMM} output.   
}
\value{
  Returns the input object invisibly.
}
\references{
  Lee, Y., Nelder, J. A. and Pawitan, Y. (2006). Generalized linear models with random effects: unified analysis via h-likelihood. Chapman & Hall: London.
}
\examples{
data("blackcap")
fit <- fitme(migStatus ~ 1+ Matern(1|longitude+latitude),data=blackcap,
             fixed=list(lambda=1,nu=1,rho=1))
plot(fit)

# compare Q-Q plots for count data
set.seed(123)
dat <- data.frame(
  x=(x <- seq(1,2,1/100)),
  y=rpois(101, lambda=x)
)
fit <- fitme(y ~ x, family=poisson(), data=dat)
plot(fit, res_type="RQR", which="mean") 
plot(fit, which="mean")
}
\keyword{hplot}
