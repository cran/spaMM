\name{Loaloa}
\alias{Loaloa}
\docType{data}
\title{
  Loa loa prevalence in North Cameroon, 1991-2001
}
\description{
This data set describes prevalence of infection by the nematode \emph{Loa loa} in North Cameroon, 1991-2001.
This is a superset of the data discussed by Diggle and Ribeiro (2007) and Diggle et al. (2007).
The study investigated the relationship between altitude, vegetation indices, and prevalence of  the parasite. 
}
\usage{data("Loaloa")}
\format{
  The data frame includes 197 observations on the following variables:
  \describe{
    \item{latitude}{latitude, in degrees.}
    \item{longitude}{longitude, in degrees.}
    \item{ntot}{sample size per location}
    \item{npos}{number of infected individuals per location}
    \item{maxNDVI}{maximum normalised-difference vegetation index (NDVI) from repeated satellite scans}
    \item{seNDVI}{standard error of NDVI}
    \item{elev1}{altitude, in m.}
    \item{elev2,elev3,elev4}{Additional altitude variables derived from the previous one, provided for convenience: 
           respectively, positive values of altitude-650, positive values of altitude-1000, and positive values of altitude-1300}
    \item{maxNDVI1}{a copy of maxNDVI modified as \code{maxNDVI1[maxNDVI1>0.8] <- 0.8}}
  }
}
\source{
The data were last retrieved on March 1, 2013 from P.J. Ribeiro's web resources
at \url{www.leg.ufpr.br/doku.php/pessoais:paulojus:mbgbook:datasets}.
}
\references{
Diggle, P., and Ribeiro, P. 2007. Model-based geostatistics, Springer series
in statistics, Springer, New York.

Diggle, P. J., Thomson, M. C., Christensen, O. F., Rowlingson, B., Obsomer,
V., Gardon, J., Wanji, S., Takougang, I., Enyong, P., Kamgno, J., Remme,
J. H., Boussinesq, M., and Molyneux, D. H. 2007. Spatial modelling and
the prediction of Loa loa risk: decision making under uncertainty, Ann.
Trop. Med. Parasitol. 101, 499-509.
}
\examples{

data("Loaloa")

### Variations on the model fit by Diggle et al. 
###    on a subset of the Loaloa data
### In each case this shows the slight differences in syntax,
###    and the difference in 'typical' computation times, 
###    when fit using corrHLfit() or fitme().

if (spaMM.getOption("example_maxtime")>4.2) {
  corrHLfit(cbind(npos,ntot-npos)~elev1+elev2+elev3+elev4+maxNDVI1+seNDVI
                   +Matern(1|longitude+latitude),HLmethod="HL(0,1)",
                 data=Loaloa,family=binomial(),ranFix=list(nu=0.5)) 
}
if (spaMM.getOption("example_maxtime")>1.8) {
  fitme(cbind(npos,ntot-npos)~elev1+elev2+elev3+elev4+maxNDVI1+seNDVI
                   +Matern(1|longitude+latitude),method="HL(0,1)",
                 data=Loaloa,family=binomial(),fixed=list(nu=0.5)) 
}

if (spaMM.getOption("example_maxtime")>7.4) {
  corrHLfit(cbind(npos,ntot-npos)~elev1+elev2+elev3+elev4+maxNDVI1+seNDVI
            +Matern(1|longitude+latitude),
              data=Loaloa,family=binomial(),ranFix=list(nu=0.5))  
}
if (spaMM.getOption("example_maxtime")>3.9) {
  fitme(cbind(npos,ntot-npos)~elev1+elev2+elev3+elev4+maxNDVI1+seNDVI
            +Matern(1|longitude+latitude),
              data=Loaloa,family=binomial(),fixed=list(nu=0.5),method="REML")
}

## Diggle and Ribeiro (2007) assumed (in this package notation) Nugget=2/7:
if (spaMM.getOption("example_maxtime")>7.2) {
  corrHLfit(cbind(npos,ntot-npos)~elev1+elev2+elev3+elev4+maxNDVI1+seNDVI
           +Matern(1|longitude+latitude),
             data=Loaloa,family=binomial(),ranFix=list(nu=0.5,Nugget=2/7))  
}
if (spaMM.getOption("example_maxtime")>2.8) {
  fitme(cbind(npos,ntot-npos)~elev1+elev2+elev3+elev4+maxNDVI1+seNDVI
           +Matern(1|longitude+latitude),method="REML",
             data=Loaloa,family=binomial(),fixed=list(nu=0.5,Nugget=2/7))  
}

## with nugget estimation:
if (spaMM.getOption("example_maxtime")>18.1) {
  corrHLfit(cbind(npos,ntot-npos)~elev1+elev2+elev3+elev4+maxNDVI1+seNDVI
           +Matern(1|longitude+latitude),
             data=Loaloa,family=binomial(),
             init.corrHLfit=list(Nugget=0.1),ranFix=list(nu=0.5))  
}
if (spaMM.getOption("example_maxtime")>7.9) {
  fitme(cbind(npos,ntot-npos)~elev1+elev2+elev3+elev4+maxNDVI1+seNDVI
           +Matern(1|longitude+latitude),
             data=Loaloa,family=binomial(),method="REML",
             init=list(Nugget=0.1),fixed=list(nu=0.5))  
}

}
\keyword{datasets}
