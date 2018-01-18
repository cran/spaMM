\name{scotlip}
\alias{scotlip}
\alias{Nmatrix}
\docType{data}
\title{Lip cancer in Scotland 1975 - 1980}
\description{
  This data set provides counts of lip cancer diagnoses made in Scottish districts from 1975 to 1980, 
  and additional information relative to these data from Clayton and Kaldor (1987) and Breslow and Clayton (1993).
  The data set contains (for each district) counts of disease events and estimates of the fraction of the population involved 
  in outdoor industry (agriculture, fishing, and forestry) which exposes it to sunlight. 

  \code{data("scotlip")} actually loads a data frame, 
  \code{scotlip}, and an adjacency matrix, \code{Nmatrix}, between 56 Scottish districts, 
  as given by Clayton and Kaldor (1987, Table 1). 
  
}
\usage{data("scotlip")}
\format{
  The data frame includes 56 observations on the following 7 variables:
  \describe{
    \item{gridcode}{alternative district identifier.}
    \item{id}{numeric district identifier (1 to 56).}
    \item{district}{district name.}
    \item{cases}{number of lip cancer cases diagnosed 1975 - 1980.}
    \item{population}{total person years at risk 1975 - 1980.}
    \item{prop.ag}{percent of the population engaged in outdoor industry.}
    \item{expec}{offsets considered by Breslow and Clayton (1993, Table 6, 'Exp' variable)}
  }
  The rows are ordered according to \code{gridcode}, so that they match the rows of \code{Nmatrix}.
}
\references{
Clayton D, Kaldor J (1987). Empirical Bayes estimates of age-standardized relative risks for use in disease mapping. Biometrics, 43: 671 - 681. 

Breslow, NE, Clayton, DG. (1993). Approximate Inference in Generalized Linear Mixed Models.
Journal of the American Statistical Association: 88 9-25.
}
\examples{
## see 'help(autoregressive)' for several exampels involving 'scotlip'.
}
\keyword{datasets}
