\name{mapMM}
\alias{mapMM}
\alias{filled.mapMM} 

\title{
  Colorful plots of predicted responses in two-dimensional space. 
}
\description{
  These functions provide either a map of predicted response in analyzed locations, or a predicted surface. The former is a straightforward representation of the analysis of the data, while the latter cope with the fact that all predictor variables may not be known in all locations on a fine spatial grid, but may involve questionable choices as a result (see \code{map.formula} argument). 
}
\usage{
mapMM(fitobject,coordinates=NULL,xrange=NULL,yrange=NULL,
      margin=1/20,add.map= FALSE, nlevels = 20, 
      color.palette = spaMM.colors,map.asp=NULL,
      col = color.palette(length(levels) - 1), plot.title, 
      plot.axes, decorations, key.title, key.axes, xaxs = "i", 
      yaxs = "i", las = 1, axes = TRUE, frame.plot = axes, ...) 
filled.mapMM(fitobject,coordinates,xrange=NULL,yrange=NULL,
             margin=1/20,map.formula,phi=1e-05,gridSteps=41,
             add.points=quote(points(pred[,coordinates],cex=1,lwd=2)),
             add.map=FALSE,axes = TRUE, plot.axes,map.asp=NULL,...)
}
\arguments{
  \item{fitobject}{
  The return object of a corrHLfit call.
%~ adjMatrix valid for mapMM ?
}
  \item{coordinates}{
  The geographical coordinates. By default they are deduced from the model formula. For example if this formula is \code{resp ~ 1 + Matern(1| x + y )} the default coordinates are c("x","y"). If this formula is \code{resp ~ 1 + Matern(1| x + y + z )}, the user must choose two of the three coordinates.
}
 \item{xrange}{
   The x range of the plot (a vector of length 2); by default defined to cover all analyzed points.  
 }
 \item{yrange}{
   The y range of the plot (a vector of length 2); by default defined to cover all analyzed points.  
 }
 \item{margin}{
   This controls how far (in relative terms) the plot extends beyond the x and y ranges of the analyzed points, and is overriden by explicit \code{xrange} and \code{yrange} arguments.  
 }
\item{map.formula}{Plotting a filled contour generally requires prediction in non-oberved locations, where predictor variables used in the original data analysis may be missing. In that case, the original model formula cannot be used and an alternative \code{map.formula} must be used to interpolate (not smooth) the predicted values in observed locations (these predictions still resulting from the original analysis based on predictor variables). As a result (1) \code{filled.mapMM} will be slower than a mere plotting function, since it involves the analysis of spatial data; (2) the results may have little useful meaning if the effects of the original predictor variables is not correctly represented by this interpolation step. For example, it may involve biases analogous to predicting temperature in non-oberved locations while ignoring effect of variation in altitude in such locations.     
}
\item{phi}{  
  This controls the phi valu assumed in the interpolation step. Ideally \code{phi} would be zero, but problems with numerically singular matrices may arise when \code{phi} is too small.
}
\item{gridSteps}{The number of levels of the grid of x and y values}
\item{add.points}{Either a boolean or an explicit expression, enclosed in \code{quote} (the default value illustrates the latter syntax). This controls whether and how analyzed locations are represented on the map.}
\item{add.map}{
  Either a boolean or an explicit expression, enclosed in \code{quote} (see Examples).
  If \code{TRUE}, the \code{map} function from the \code{maps} package (which much therefore the loaded) is used to add a map from its default \code{world} database. \code{xrange} and \code{yrange} are used to select the area, so it is most convenient if the \code{coordinates} are longitude and latitude (in this order and in standard units). An explicit expression can also be used for further control.  
}
\item{levels}{  
a set of levels which are used to partition the range of z. Must be strictly increasing (and finite). Areas with z values between consecutive levels are painted with the same color.
}
\item{nlevels}{	
if levels is not specified, the range of z, values is divided into approximately this many levels.
}
\item{color.palette}{	
a color palette function to be used to assign colors in the plot.
}
\item{map.asp}{  
the y/x aspect ratio of the 2D plot area (not of the full figure including the scale). Default is (plotted y range)/(plotted x range) (i.e., scales for x are identical). 
}
\item{col}{	
an explicit set of colors to be used in the plot. This argument overrides any palette function specification. There should be one less color than levels
}
\item{plot.title}{	
statements which add titles to the main plot. If provided, further \code{...} arguments are ignored (see Details).
}
\item{plot.axes}{	
statements which draw axes (and a box) on the main plot. Default axes are drawn by default when this argument is missing, given \code{axes = TRUE}.
}
\item{decorations}{  
Additional graphic statements (\code{points}, \code{polygon}, etc.).
}
\item{key.title}{	
statements which add titles for the plot key.
}
\item{key.axes}{	
statements which draw axes on the plot key. 
}
\item{xaxs}{	
the x axis style. The default is to use internal labeling.
}
\item{yaxs}{	
the y axis style. The default is to use internal labeling.
}
\item{las}{	
the style of labeling to be used. The default is to use horizontal labeling.
}
\item{axes, frame.plot}{	
logicals indicating if axes and a box should be drawn, as in plot.default.
}
\item{\dots}{
  further arguments passed to or from other methods. For \code{mapMM}, currently only additional graphical parameters passed to \code{title()} (see Details).
  For \code{filled.mapMM}, these parameters are those that can be passed to \code{\link{spaMM.filled.contour}}. 
}
}
\value{
  No return value. Plots are produced as side-effects. 
}
\details{
If you have values for all predictor variables in all locations of a fine spatial grid, \code{filled.mapMM} may not be a good choice.
Rather, use \code{predict(<fitobject>,newX= <all predictor variables >)} to generate all predictions, and then either 
\code{spaMM.filled.contour} or some other raster functions.

These functions handle some of their arguments as \code{\link{filled.contour}} does. For \code{mapMM} in particular, this means that either
\code{plot.title} is missing, or \code{...} is ignored. Thus, one can provide an optional \code{xlab} either in the \code{plot.title} argument, or in the \code{...} \code{plot.title} is missing. \code{filled.mapMM} calls \code{spaMM.filled.contour} which behaves identically, so the \code{...} argument of \code{filled.mapMM} should contain either a \code{plot.title} or further arguments. 

A side effet is that \code{filled.mapMM}, like \code{filled.contour}, does not does not provide axis labels (\code{xlab} and \code{ylab}) by default. 
}

\examples{
data(blackcap)
bfit <- corrHLfit(migStatus ~ means+ Matern(1|longitude+latitude),data=blackcap,
                  HLmethod="ML",
                  ranFix=list(lambda=0.5537,phi=1.376e-05,rho=0.0544740,nu=0.6286311))
mapMM(bfit,color.palette = function(n){spaMM.colors(n,redshift=1/2)})

\donttest{
## showing add.map
filled.mapMM(bfit,add.map=TRUE,plot.axes={axis(1);axis(2)},
             plot.title=title(main="Inferred migration propensity of blackcaps",
                               xlab="longitude",ylab="latitude"))
## filled.mapMM takes a bit longer
filled.mapMM(bfit,nlevels=30,plot.axes={axis(1);axis(2)},
              plot.title=title(main="Inferred migration propensity of blackcaps",
                               xlab="longitude",ylab="latitude"))
}


## showing more involved use of add.map
data(Loaloa)  
lfit <- corrHLfit(cbind(npos,ntot-npos)~elev1+elev2+elev3+elev4+maxNDVI1+seNDVI
                    +Matern(1|longitude+latitude),HLmethod="HL(0,1)",data=Loaloa,
                    family=binomial(),ranFix=list(nu=0.5,rho=2.255197,lambda=1.075))   

if(require(maps)) {
  mapMM(lfit,add.map=quote(map(,xlim=xrange,ylim=yrange,add=TRUE)))
}


\donttest{
## longer computation requiring interpolation of 197 points 
filled.mapMM(lfit,add.map=TRUE,plot.axes={axis(1);axis(2)},
             add.points=quote(points(pred[,coordinates],pch=15,cex=0.3)),
             plot.title=title(main="Inferred prevalence, North Cameroon",
                                xlab="longitude",ylab="latitude"))
}

}
\keyword{hplot}