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
\usage{data(Loaloa)}
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
    \item{maxNDVI1}{\code{maxNDVI1[maxNDVI1>0.8] <- 0.8}}
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
data(Loaloa)
## Variations on the model fit by Diggle et al. on a subset of these data

## With HLmethod=HL(0,1), this takes ~35s on a recent laptop
\dontrun{
corrHLfit(cbind(npos,ntot-npos)~elev1+elev2+elev3+elev4+maxNDVI1+seNDVI
                   +Matern(1|longitude+latitude),HLmethod="HL(0,1)",
                 data=Loaloa,family=binomial(),ranFix=list(nu=0.5)) 
## nu=0.5 is the exponential correlation     
}

## the following fits take ~ 1 -- 5 minutes each.

\dontrun{
corrHLfit(cbind(npos,ntot-npos)~elev1+elev2+elev3+elev4+maxNDVI1+seNDVI
            +Matern(1|longitude+latitude),
              data=Loaloa,family=binomial(),ranFix=list(nu=0.5)) 
}

## Diggle and Ribeiro (2007) assumed (in this package notation) Nugget=2/7:
\dontrun{
corrHLfit(cbind(npos,ntot-npos)~elev1+elev2+elev3+elev4+maxNDVI1+seNDVI
           +Matern(1|longitude+latitude),
             data=Loaloa,family=binomial(),ranFix=list(nu=0.5,Nugget=2/7))
}

## with nugget estimation:
\dontrun{
corrHLfit(cbind(npos,ntot-npos)~elev1+elev2+elev3+elev4+maxNDVI1+seNDVI
           +Matern(1|longitude+latitude),
             data=Loaloa,family=binomial(),
             init.corrHLfit=list(Nugget=0.1),ranFix=list(nu=0.5)) 
}

}
\keyword{datasets}
