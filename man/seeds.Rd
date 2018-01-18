\name{seeds}
\alias{seeds}
\docType{data}
\title{
Seed germination data
}
\description{
A classic toy data set, \dQuote{from research conducted by microbiologist Dr P. Whitney of Surrey University. A batch of tiny seeds is brushed onto a plate covered with a certain extract at a given dilution. The numbers of germinated and ungerminated seeds are subsequently counted} (Crowder, 1978). Two seed types and two extracts are here considered in a 2x2 factorial design.
}
\usage{data("seeds")}
\format{
  The data frame includes 21 observations on the following variables:
  \describe{
    \item{plate}{Factor for replication;}
    \item{seed}{Seed type, a factor with two levels O73 and O75;}
    \item{extract}{Root extract, a factor with two levels Bean and Cucumber;}
    \item{r}{Number of seeds that germinated;}
    \item{n}{Total number of seeds tested}
  }
}
\source{
Crowder (1978), Table 3.
}
\references{
Crowder, M.J., 1978. Beta-binomial anova for proportions. Appl. Statist., 27, 34-37. 

Y. Lee and J. A. Nelder. 1996. Hierarchical generalized linear models (with discussion). J. R. Statist. Soc. B, 58: 619-678. 
}
\examples{
data("seeds") 
## An extended quasi-likelihood (EQL) fit as considered by Lee and Nelder (1996):
HLfit(cbind(r,n-r)~seed*extract+(1|plate),family=binomial(),
  rand.family=Beta(),
  HLmethod="HL(0,0)",
  data=seeds)
}
\keyword{datasets}
