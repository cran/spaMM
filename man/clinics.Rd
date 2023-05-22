\name{clinics}
\alias{clinics}
\docType{data}
\title{
Toy dataset for binomial response
}
\description{
A small data set used by Booth & Hobert (1998). 
}
\usage{data("clinics")}
\format{
  A data frame with 16 observations on the following 4 variables.
  \describe{
    \item{\code{npos}}{a numeric vector}
    \item{\code{nneg}}{a numeric vector}
    \item{\code{treatment}}{a numeric vector}
    \item{\code{clinic}}{a numeric vector}
  }
}
\references{
  Booth, J.G., Hobert, J.P. (1998) Standard errors of prediction in generalized linear mixed models. J. Am. Stat. Assoc. 93: 262-272. 
}
\examples{
data(clinics)
\dontrun{
# The dataset was built as follows
npos <- c(11,16,14,2,6,1,1,4,10,22,7,1,0,0,1,6)
ntot <- c(36,20,19,16,17,11,5,6,37,32,19,17,12,10,9,7)
treatment <- c(rep(1,8),rep(0,8))
clinic <-c(seq(8),seq(8))
clinics <- data.frame(npos=npos,nneg=ntot-npos,treatment=treatment,clinic=clinic)
}
}
\keyword{datasets}
