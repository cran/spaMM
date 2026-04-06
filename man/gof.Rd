\name{gof}
\alias{gof}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Goodness of fit test
}
\description{
Performs a test of goodness of fit (by default: see the \code{testfn} argument for possible alternatives). The only method implemented so far is based on randomized quantile residuals (Dunn & Smyth, 1996).
}
\usage{
gof(object, method = "RQR", plot.=FALSE, testfn=stats::shapiro.test, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
A fit object inheriting from class \code{"HLfit"}, as returned by \pkg{spaMM}'s main fitting functions.
}
  \item{method}{
Character string; a method of test of goodness of fit.
}
\item{plot.}{Boolean: whether to produce a plot (effectively \code{plot(object, res_type="RQR", which="mean")} for RQRs as function of predicted values).}
\item{testfn}{ A function with the following interface: it takes the vector of values whose normality is assessed as first argument; and it returns a \code{list} (to which the \code{RQR} element will be added, see Value). There is no particular constraint on the contents of the returned list. For example, \code{\link[stats]{qqnorm}} is a possible \code{testfn} although it does not return a test.}
\item{\dots}{
    Further arguments passed to \code{testfn}. 
}
}
\value{
Returns a list made up of the return value of a call to the function provided by the \code{testfn} argument, with added element \code{RQR}, itself a list including the randomized residuals.
The default \code{testfn}, \code{\link[stats]{shapiro.test}} (except for datasets of more than 5000 fitted response values), returns \dQuote{a list with class \code{"htest"}} which is therefore the class of the default result. 
}
\references{
Dunn, K. P., and Smyth, G. K. (1996). Randomized quantile residuals. Journal of Computational and Graphical Statistics 5, 1-10.
}
\examples{
data("Orthodont", package = "nlme")
sp1 <- fitme(distance ~ age * Sex + (1 | Subject), data = Orthodont)
gof(sp1) #  classic toy example, but poor fit.                
gof(sp1, testfn=ks.test, y="pnorm")                
}
\keyword{htest}