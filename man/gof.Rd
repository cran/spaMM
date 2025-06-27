\name{gof}
\alias{gof}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Goodness of fit test
}
\description{
Performs a test of goodness of fit. The only method implemented so far is based on randomized quantile residuals (Dunn & Smyth, 1996).
}
\usage{
gof(object, method = "RQR", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
A fit object inheriting from class \code{"HLfit"}, as returned by \pkg{spaMM}'s main fitting functions.
}
  \item{method}{
Character string; a method of test of goodness of fit.
}
  \item{\dots}{
    Further arguments possibly passed to internal methods. 
}
}
\value{
Returns the value of a call to \code{\link[stats]{shapiro.test}} (\dQuote{a list with class \code{"htest"}}) with added element \code{RQR}, itself a list including the randomized residuals.
}
\references{
Dunn, K. P., and Smyth, G. K. (1996). Randomized quantile residuals. Journal of Computational and Graphical Statistics 5, 1-10.
}
\examples{
data("Orthodont", package = "nlme")
sp1 <- fitme(distance ~ age * Sex + (1 | Subject), data = Orthodont)
gof(sp1) #  classic toy example, but poor fit.                
}
\keyword{htest}