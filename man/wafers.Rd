\name{wafers}
\alias{wafers}
\docType{data}
\title{Data from a resistivity experiment for semiconductor materials.}
\description{
  This data set was reported and analyzed by Robinson et al. (2006) and reanalyzed by Lee et al. (2011). The data 
  \dQuote{deal with wafers in a single etching process in semiconductor manufacturing.
Wafers vary through time since there are some variables that are not perfectly controllable in the etching process. For this
reason, wafers produced on any given day (batch) may be different from those produced on another day (batch). To measure
variation over batch, wafers are tested by choosing several days at random. In this data, resistivity is the response of interest.
There are three variables, gas flow rate (x1), temperature (x2), and pressure (x3) and one random effect (batch or day).}  (Lee et al 2011).
}
\source{
This data set was manually pasted from Table 3 of Lee et al. (2011). Transcription errors may have occurred. 
}
\usage{data("wafers")}
\format{
  The data frame includes 198 observations on the following variables:
  \describe{
    \item{y}{resistivity.}
    \item{batch}{batch, indeed.}
    \item{X1}{gas flow rate.}
    \item{X2}{temperature.}
    \item{X3}{pressure.}
  }
}
\references{
Robinson TJ, Wulff SS, Montgomery DC, Khuri AI. 2006. Robust parameter design using generalized linear mixed models. Journal of Quality
Technology 38: 38--65.

Lee, Y., Nelder, J.A., and Park, H. 2011. HGLMs for quality improvement. 
     Applied Stochastic Models in Business and Industry 27, 315-328.
}
\examples{
## see examples in the main Documentation page for the package.
}
\keyword{datasets}
