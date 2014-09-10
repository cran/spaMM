\name{seaMask}
\alias{seaMask}
\alias{landMask}
\docType{data}
\title{
  Masks of seas or lands
}
\description{
  These convenient masks can be added to maps of (parts of) the world to mask map information for these areas. 
}
\usage{
data(seaMask)
data(landMask)
}
\format{
\code{seaMask} and \code{landMask} are data frames with two variables, \code{x} and \code{y} for longitude and latitude. 
Its contents are suitable for use with \code{\link{polypath}}: they define different polygones, each separated by a row of \code{NA}s.
}
\details{
\code{polypath} requires polygons, while \code{map(interior=FALSE,plot=FALSE)} returns small segments. \code{landMask} is the result of reconnecting the segments into full coastlines of all land blocks. 
}
\examples{
data(seaMask)
## plot of predictions of behaviour for a land bird:
if (require(maps)){
 data(blackcap)
 bfit <- corrHLfit(migStatus ~ means+ Matern(1|longitude+latitude),data=blackcap,
                  HLmethod="ML",
                  ranFix=list(lambda=0.5537,phi=1.376e-05,rho=0.0544740,nu=0.6286311))
 ## We add small masks to the points on small islands to see the predictions there
 ll <- blackcap[,c("longitude","latitude")]
 pointmask <- function(xy,r=1,npts=12) {
  theta <- 2*pi/npts *seq(npts)
  hexas <- lapply(seq(nrow(xy)),function(li){
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
             add.points=quote(points(pred[,coordinates],cex=1,pch="+")),
             plot.axes={axis(1);axis(2);
                        polypath(rbind(seaMask,pmasks),border=FALSE,
                                 col="grey", rule="evenodd")
             })
}
}
\keyword{datasets}
