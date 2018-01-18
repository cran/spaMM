\name{arabidopsis}
\alias{arabidopsis}
\docType{data}
\title{
  Arabidopsis genetic and climatic data
}
\description{
  For 948 \dQuote{accessions} from European Arabidopsis thaliana populations, this data set merges the genotypic information at 
  four single nucleotide polmorphisms (SNP) putatively involved in adaptation to climate (Fournier-Level et al, 2011, Table 1), 
  with 13 climatic variables from Hancock et al. (2011). 
}
\usage{data("arabidopsis")}
\format{
  The data frame includes 948 observations on the following variables:
  \describe{
    \item{pos1046738, pos5510910, pos6235221, pos8132698}{Genotypes at four SNP loci}
    \item{LAT}{latitude}
    \item{LONG}{longitude}
    \item{seasonal, tempWarmest, tempColdest, preciWettest,
       preciDriest, preciCV, PAR_SPRING,}{} 
    \item{growingL, conseqCold, conseqFrFree, RelHumidSp, dayLSp, aridity}{Thirteen climatic variables. 
    See Hancock et al. (2011) for details about these variables.}
  }
}
\source{
The data were retrieved from \code{http://bergelson.uchicago.edu/regmap-data/climate-genome-scan} on 22 February 2013 (they may no longer be available from there).
}
\details{
  The  response is binary so \code{HLmethod="PQL/L"} seems warranted (see Rousset and Ferdy, 2014).  
}
\references{
 Fournier-Level A, Korte A., Cooper M. D., Nordborg M., Schmitt J., Wilczek AM (2011). A map of local adaptation in Arabidopsis thaliana. Science 334: 86-89.

 Hancock, A. M., Brachi, B., Faure, N., Horton, M. W., Jarymowycz, L. B., Sperone, F. G., Toomajian, C., Roux, F., and Bergelson, J. 2011. 
 Adaptation to climate across the Arabidopsis thaliana genome, Science 334: 83-86.
 
  Rousset F., Ferdy, J.-B. (2014) Testing environmental and genetic effects in the presence of spatial autocorrelation. Ecography, 37: 781-790.
    \url{http://dx.doi.org/10.1111/ecog.00566}
 
}
\examples{
data("arabidopsis")
if (spaMM.getOption("example_maxtime")>8) {
  HLCor(cbind(pos1046738,1-pos1046738)~seasonal+Matern(1|LAT+LONG),
        ranPars=list(rho=0.1192779,nu=0.2369892,lambda=8.599),
        family=binomial(),HLmethod="PQL/L",data=arabidopsis)
}
## The above ranPars are deduced from the following fit:
if (spaMM.getOption("example_maxtime")>46) {
  SNPfit <- fitme(cbind(pos1046738,1-pos1046738)~seasonal+Matern(1|LAT+LONG),
              verbose=c(TRACE=TRUE), 
              family=binomial(),method="PQL/L",data=arabidopsis)
  summary(SNPfit) # p_v=-125.0392
}
}
\keyword{datasets}
