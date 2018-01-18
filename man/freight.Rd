\name{freight}
\alias{freight}
\docType{data}
\alias{freight}
\title{Freight dataset}
\description{
  A set of data on airfreight breakage. Data are given on 10 air shipments, each carrying
1000 ampules of some substance. For each shipment, the number of ampules found broken upon arrival, and the number of times the shipments were transferred from one aircraft to another, are recorded.
}
\usage{
data("freight")
}
\format{
 The data frame includes 10 observations on the following variables:
  \describe{
    \item{broken}{number of ampules found broken upon arrival.}
    \item{transfers}{number of times the shipments were transferred from one aircraft to another.}
    \item{id}{Shipment identifier.}
  }
}
\source{
The data set is reported by Kutner et al. (2003) and used by Sellers and Shmueli (2010) to illustrate COMPoisson analyses.
}
\references{
	Kutner MH, Nachtsheim CJ, Neter J, Li W (2005, p. 35). Applied Linear Regression Models, Fourth Edition. McGraw-Hill.
	
	Sellers KF, Shmueli G (2010) A Flexible Regression Model for Count Data. Ann. Appl. Stat. 4: 943–961
}
\keyword{datasets}
\examples{
## see ?COMPoisson for examples
}