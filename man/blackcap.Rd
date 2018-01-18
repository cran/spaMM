\name{blackcap}
\alias{blackcap}
\docType{data}
\title{
Genetic polymorphism in relation to migration in the blackcap
}
\description{
This data set is extracted from a study of genetic polymorphisms potentially associated to migration behaviour in the blackcap
(Sylvia atricapilla). Across different populations in Europe and Africa, the average migration behaviour 
was found to correlate with average allele size (dependent on the number of repeats of a small DNA motif) 
at the locus ADCYAP1, encoding a neuropeptide.
This data set is quite small and ill-suited for separating random-effect variance from residual variance. The likelihood surface for the Matérn model actually has local maxima. 
}
\usage{data("blackcap")}
\format{
  The data frame includes 14 observations on the following variables:
  \describe{
    \item{latitude}{latitude, indeed.}
    \item{longitude}{longitude, indeed.}
    \item{migStatus}{migration status as determined by Mueller et al, from 0 (resident populations) to 2.5 (long-distance migratory populations)}
    \item{means}{Mean allele sizes in each population}
    \item{pos}{Numerical index for the populations}
  }
}
\details{
Migration status was coded as : pure resident populations as '0', resident populations with 
some migratory restlessness as '0.5', partial migratory populations as '1', completely migratory populations migrating
short-distances as '1.5', intermediate-distance migratory populations as '2' and distinct long-distance migratory populations
as '2.5'.
}
\source{
Data from Mueller et al. (2011), including supplementary material retrieved from
\url{http://rspb.royalsocietypublishing.org/content/suppl/2011/02/11/rspb.2010.2567.DC1}.
}
\references{
Mueller, J. C., Pulido, F., and Kempenaers, B. 2011. Identification of a
gene associated with avian migratory behaviour, Proc. Roy. Soc. (Lond.)
B 278, 2848-2856.
}
\examples{
## see 'corrHLfit' and 'fixedLRT' for examples involving these data
}
\keyword{datasets}
