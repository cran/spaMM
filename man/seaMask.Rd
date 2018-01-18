\name{seaMask}
\alias{seaMask}
\alias{oceanmask}
\alias{landMask}
\alias{worldcountries}
\docType{data}
\title{
  Masks of seas or lands
}
\description{
  These convenient masks can be added to maps of (parts of) the world to mask map information for these areas. 
}
\usage{
data("seaMask")
data("landMask")
data("worldcountries")
data("oceanmask")
}
\format{
\code{seaMask} and \code{landMask} are data frames with two variables, \code{x} and \code{y} for longitude and latitude. 
Its contents are suitable for use with \code{\link{polypath}}: they define different polygones, each separated by a row of \code{NA}s.

\code{worldcountries} and \code{oceanmask} are \code{SpatialPolygonsDataFrame} objects.
% hence require 'sp', but no (runnable) example require sp.
% R CMD check may complain if 'sp' is not available locally for a check (but it will be available for checks on CRAN). 
}
\details{
A land mask can be produced out of \code{worldcountries} by filling the countries (i.e. \code{fill="black"} in the code for \code{country.layer} in the Examples in \url{http://kimura.univ-montp2.fr/~rousset/spaMM/example_raster.html}.

\code{worldcountries} and \code{oceanmask} were created from public domain shapefiles downloaded from \url{www.naturalearth.com} on 2015/10/21. These are suitable for plots involving geographical projections not available through \code{map}, and more generally for raster plots. Only the lowest-resolution data are included in spaMM, to minimize the size of the package archive, but higher-resolution files are available on \url{www.naturalearth.com}, from where they can be loaded as shown in the examples. \code{worldcountries} had to be edited for non-ASCII characters before inclusion in spaMM:
\code{worldcountries@data$formal_fr} was removed and the \dQuote{Côte d'Ivoire} level of some factor variables was renamed.  

\code{seaMask} and \code{landMask} were created from the world map in the maps package. 
\code{polypath} requires polygons, while \code{map(interior=FALSE,plot=FALSE)} returns small segments. \code{landMask} is the result of reconnecting the segments into full coastlines of all land blocks. 


}
\seealso{\url{http://kimura.univ-montp2.fr/~rousset/spaMM/example_raster.html} for uses of \code{worldcountries} and \code{oceanmask}} 
\examples{
if (spaMM.getOption("example_maxtime")>1.1) {

data("seaMask")
## plot of predictions of behaviour for a land bird:
if (require(maps)){
 data("blackcap")
 bfit <- corrHLfit(migStatus ~ means+ Matern(1|longitude+latitude),data=blackcap,
                  HLmethod="ML",
                  ranFix=list(lambda=0.5537,phi=1.376e-05,rho=0.0544740,nu=0.6286311))
 ## We add small masks to the points on small islands to see the predictions there
 ll <- blackcap[,c("longitude","latitude")]
 pointmask <- function(xy,r=1,npts=12) {
  theta <- 2*pi/npts *seq(npts)
  hexas <- lapply(seq(nrow(xy)), function(li){
    p <- as.numeric(xy[li,])
    hexa <- cbind(x=p[1]+r*cos(theta),y=p[2]+r*sin(theta))
    rbind(rep(NA,2),hexa) ## inital NA before each polygon
  })
  do.call(rbind,hexas)
 }
 pmasks <- pointmask(ll[c(2,4,5,6,7),],r=0.8) ## small islands only
 filled.mapMM(bfit,add.map=TRUE,
             plot.title=title(main="Inferred migration propensity of blackcaps",
                               xlab="longitude",ylab="latitude"),
             decorations=quote(points(pred[,coordinates],cex=1,pch="+")),
             plot.axes=quote({axis(1);axis(2);
                        polypath(rbind(seaMask,pmasks),border=FALSE,
                                 col="grey", rule="evenodd")
             }))
}
}

\dontrun{
# All shape files can be found here: http://www.naturalearthdata.com/downloads/
# Once downloaded, they can be loaded into R by
if (require("rgdal", quietly = TRUE)) {
  worldcountries <- readOGR("ne_110m_admin_0_countries_lakes.shp",
                           layer="ne_110m_admin_0_countries_lakes")
}
}
}
\keyword{datasets}
