\name{mapMM}
\alias{spaMMplot2D}
\alias{mapMM}
\alias{filled.mapMM} 
\alias{map_ranef} 

\title{
  Colorful plots of predictions in two-dimensional space. 
}
\description{
  These functions provide either a map of predicted response in analyzed locations, or a predicted surface. \code{mapMM} is a straightforward representation of the analysis of the data, while \code{filled.mapMM} uses interpolation to cope with the fact that all predictor variables may not be known in all locations on a fine spatial grid.  \code{map_ranef} maps a single spatial random effect. These three functions takes an \code{HLfit} object as input. \code{mapMM} calls \code{spaMMplot2D}, which is similar but takes a more conventional (x,y,z) input. 
  
  Using \code{filled.mapMM} may involve questionable choices. Plotting a filled contour generally requires prediction in non-observed locations, where predictor variables used in the original data analysis may be missing. In that case, the original model formula cannot be used and an alternative model (controlled by the \code{map.formula} argument) must be used to interpolate (not smooth) the predicted values in observed locations (these predictions still resulting from the original analysis based on predictor variables). \code{filled.mapMM} always performs such interpolation (it does not allow one to provide values for the predictor variables). As a result (1) \code{filled.mapMM} will be slower than a mere plotting function, since it involves the analysis of spatial data; (2) the results may have little useful meaning if the effect of the original predictor variables is not correctly represented by this interpolation step. For example, prediction by interpolation may be biased in a way analogous to prediction of temperature in non-observed locations while ignoring effect of variation in altitude in such locations. Likewise, the \code{variance} argument of \code{filled.mapMM} allows one only to plot the prediction variance of its own interpolator, rather than that of the input object. 
  
  \code{map_ranef} is free of the limitations of \code{filled.mapMM}. 
}
\usage{
spaMMplot2D(x, y, z, xrange=range(x, finite = TRUE),
  yrange=range(y, finite = TRUE), margin=1/20, add.map= FALSE, 
  nlevels = 20, color.palette = spaMM.colors, map.asp=NULL,
  col = color.palette(length(levels) - 1), plot.title=NULL, plot.axes=NULL, 
  decorations=NULL, key.title=NULL, key.axes=NULL, xaxs = "i", 
  yaxs = "i", las = 1, axes = TRUE, frame.plot = axes, ...) 

mapMM(fitobject,Ztransf=NULL,coordinates,
  add.points,decorations=NULL,plot.title=NULL,plot.axes=NULL,envir=-3, ...)

filled.mapMM(
  fitobject, Ztransf = NULL, coordinates, xrange = NULL, yrange = NULL, 
  margin = 1/20, map.formula, phi = 1e-05, gridSteps = 41, 
  decorations = quote(points(pred[, coordinates], cex = 1, lwd = 2)),
  add.map = FALSE, axes = TRUE, plot.title = NULL, plot.axes = NULL, 
  map.asp = NULL, variance = NULL, var.contour.args = list(), 
  smoothObject = NULL, return.="smoothObject", ...)

map_ranef(fitobject, re.form, Ztransf=NULL, xrange = NULL, yrange = NULL, 
  margin = 1/20, gridSteps = 41, 
  decorations = quote(points(fitobject$data[, coordinates], cex = 1, lwd = 2)), 
  add.map = FALSE, axes = TRUE, plot.title=NULL, plot.axes=NULL, 
  map.asp = NULL, ...)
}
\arguments{
  \item{fitobject}{
  The return object of a corrHLfit or fitme call.
%~ adjMatrix valid for mapMM ?
}
  \item{x,y,z}{
  Three vectors of coordinates, with \code{z} being expectedly the response.
%~ adjMatrix valid for mapMM ?
}
\item{re.form}{
  A model formula giving the single random effect term to plot, needed only if there are several spatial random effects in the fitted model. In that case, it must be formatted as \code{ . ~ <term>}, as for the \code{re.form} argument of \code{predict.HLfit}.
}
\item{Ztransf}{
  A transformation of the predicted response, given as a function whose only required argument can be a one-column matrix. The name of this argument must be 
   \code{Z} (not \code{x}), as is appropriate for use in \code{do.call(Ztransf,list(Z=Zvalues))}. 
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
\item{map.formula}{NULL, or a formula whose left-hand side is ignored. Provides the formula used for interpolation. If NULL, a default formula with the same spatial effect(s) as in the input \code{fitobject} is used.     
}
\item{phi}{  
  This controls the phi value assumed in the interpolation step. Ideally \code{phi} would be zero, but problems with numerically singular matrices may arise when \code{phi} is too small.
}
\item{gridSteps}{The number of levels of the grid of x and y values}
\item{variance}{Either NULL, or the name of a component of variance of prediction \emph{by the interpolator} to be plotted. Must name one of the components that can be returned by \code{predict.HLfit}. \code{variance="predVar"} is suitable for uncertainty in point prediction.}
\item{var.contour.args}{A list of control parameters for rendering of prediction variances. See \code{\link{contour}} for possible arguments (except \code{x}, \code{y}, \code{z} and \code{add}).}
\item{add.map}{
  Either a boolean or an explicit expression, enclosed in \code{quote} (see Examples).
  If \code{TRUE}, the \code{map} function from the \code{maps} package (which much therefore the loaded) is used to add a map from its default \code{world} database. \code{xrange} and \code{yrange} are used to select the area, so it is most convenient if the \code{coordinates} are longitude and latitude (in this order and in standard units). An explicit expression can also be used for further control.  
}
\item{levels}{  
a set of levels which are used to partition the range of z. Must be strictly increasing (and finite). Areas with z values between consecutive levels are painted with the same color.
}
\item{nlevels}{	
if \code{levels} is not specified, the range of z, values is divided into *approximately* this many levels (a call to \code{\link{pretty}} determines the actual number of levels).
}
\item{color.palette}{	
a color palette function to be used to assign colors in the plot.
}
\item{map.asp}{  
the y/x aspect ratio of the 2D plot area (not of the full figure including the scale). By default, the scales for x and y are identical unless the x and y ranges are too different. Namely, the scales are identical if (plotted y range)/(plotted x range) is 1/4 < . < 4, and map.asp is 1 otherwise. 
}
\item{col}{	
an explicit set of colors to be used in the plot. This argument overrides any palette function specification. There should be one less color than levels
}
\item{plot.title}{	
statements which add titles to the main plot. See Details for differences between functions.
}
\item{plot.axes}{	
statements which draw axes (and a box) on the main plot. See Details for differences between functions. 
}
\item{decorations}{ Either NULL or Additional graphic statements (\code{points}, \code{polygon}, etc.), enclosed in \code{quote} (the default value illustrates the latter syntax).
.
}
\item{add.points}{Obsolete, use \code{decorations} instead.}
\item{envir}{Controls the environment in which \code{plot.title}, \code{plot.axes}, and 
    \code{decorations} are evaluated. \code{mapMM} calls \code{spaMM2Dplot} from where these graphic arguments are evaluated, and the default value -3 means that they are evaluated within the environment from where \code{mapMM} was called.}
\item{key.title}{	
statements which add titles for the plot key.
}
\item{key.axes}{	
statements which draw axes on the plot key. 
}

\item{xaxs}{the x axis style. The default is to use internal labeling.}

\item{yaxs}{the y axis style. The default is to use internal labeling.}

\item{las}{the style of labeling to be used. The default is to use horizontal labeling.}

\item{axes, frame.plot}{logicals indicating if axes and a box should be drawn, as in plot.default.}

\item{smoothObject}{Either NULL, or an object inheriting from class \code{HLfit} (hence, an object on which \code{predict.HLfit} can be called), predicting the response surface in any coordinates. See Details for typical usages.}
\item{return.}{character string: see Value}

\item{\dots}{
  further arguments passed to or from other methods. For \code{mapMM}, all such arguments are passed to \code{spaMMplot2D}; for \code{spaMMplot2D}, currently only additional graphical parameters passed to \code{title()} (see Details).  For \code{filled.mapMM} and \code{map_ranef}, these parameters are those that can be passed to \code{\link{spaMM.filled.contour}}. 
}
}
\value{
  \code{filled.mapMM} by default returns invisibly the fit object predicting the interpolated response surface; however, for any non-default \code{return.} argument (\code{return.="raster"} would be recommended to ensure future back-compatibility), it will return a raster of values as a list with elements \code{x}, \code{y} and \code{z}. \code{map_ranef} returns invisibly a 3-column matrix containing the spatial coordinates, and the predicted effect \code{z} on the linear predictor scale (which is also the scale of the plot, unless a \code{Ztransf} is used). \code{mapMM} returns invisibly a list with elements \code{x}, \code{y} and \code{z}. Plots are produced as side-effects. 
}
\details{

The \code{smoothObject} argument may be used to redraw a figure faster by recycling the predictor of the response surface returned invisibly by a previous call to \code{filled.mapMM}.   

For \code{smoothObject=NULL} (the default), \code{filled.mapMM} interpolates the predicted response, with sometimes unpleasant effects. For example, if one interpolates probabilities, the result may not be within [0,1], and then (say) a logarithmic \code{Ztransf} may generate NaN values that would otherwise not occur. The \code{smoothObject} argument may be used to overcome the default behaviour, by providing an alternative predictor. 

If you have values for all predictor variables in all locations of a fine spatial grid, \code{filled.mapMM} may not be a good choice, since it will ignore that information (see \code{map.formula} argument). Rather, one should use \code{predict(<fitobject>,newdata= <all predictor variables >)} to generate all predictions, and then either 
\code{spaMM.filled.contour} or some other raster functions.

%The functions handle some of their arguments as \code{\link{filled.contour}} does. For \code{mapMM} in particular, this means that either \code{plot.title} is NULL, or \code{...} is ignored. Thus, one can provide an optional \code{xlab} either in the \code{plot.title} argument, or in the \code{...} if \code{plot.title} is missing. \code{filled.mapMM} calls \code{spaMM.filled.contour} which behaves identically, so the \code{...} argument of \code{filled.mapMM} should contain either a \code{plot.title} or further arguments. 

The different functions are (currently) inconsistent among themselves in the way they handle the \code{plot.title} and \code{plot.axes} argument:

\bold{spaMM.filled.contour} behaves like \code{graphics::filled.contour}, which (1) handles arguments which are calls such as \code{title(.)} or \code{{axis(1);axis(2)}}; (2) ignores \code{...} arguments if \code{plot.title} is missing; and (3) draws axes by default when \code{plot.axes} is missing, given \code{axes = TRUE}.

By contrast, \bold{filled.mapMM} handles arguments which are language expressions such as produced by \code{quote(.)} or \code{substitute(.)} (see Examples). 

\bold{mapMM} can handles language expressions, but also accepts at least some calls.
}
\seealso{
\code{\link{seaMask}} for masking areas in a filled map;
\url{https://gitlab.mbb.univ-montp2.fr/francois/spamm-ref/-/blob/master/vignettePlus/example_raster.html} for more elaborate plot procedures.
}
\examples{
data("blackcap")
bfit <- fitme(migStatus ~ means+ Matern(1|longitude+latitude),data=blackcap,
              fixed=list(lambda=0.5537,phi=1.376e-05,rho=0.0544740,nu=0.6286311))
mapMM(bfit,color.palette = function(n){spaMM.colors(n,redshift=1/2)},add.map=TRUE)
map_ranef(bfit) # providing argument re.form= . ~ Matern(1|longitude+latitude)

if (spaMM.getOption("example_maxtime")>1) {
  ## filled.mapMM takes a bit longer
  # showing 'add.map', 'nlevels', and contour lines for 'variance'
  filled.mapMM(bfit, nlevels=30, add.map=TRUE, plot.axes=quote({axis(1);axis(2)}),
             variance="respVar",
             plot.title=title(main="Inferred migration propensity of blackcaps",
                               xlab="longitude",ylab="latitude"))
                               
  ## Similar plots by ggplot2:
  \dontrun{
   library(rnaturalearth) # provides sea mask through 'ne_download' function
   library(ggplot2)
   library(sp)
  
   # sea mask 
   sea <- ne_download(scale = 10, type = 'ocean', category = "physical", returnclass = "sf")
   
   # Generation of data.frame for ggplot:
   rastr <- filled.mapMM(bfit, return.="raster")
   spdf <- data.frame(Long=rep(rastr$x, nc), Lat=rastr$y[gl(nr,nc)], z = as.vector(rastr$z))

   ggplot(spdf) + 
     geom_contour_filled(aes(Long,Lat,z=z), bins = 20) +
     guides(fill = "none") +
     geom_sf(data = sea, fill = "black") +
     coord_sf(ylim = range(rastr$y), xlim = range(rastr$x), expand = FALSE)
  }
                               
}

if (spaMM.getOption("example_maxtime")>3) {
 data("Loaloa")  
 lfit <- fitme(cbind(npos,ntot-npos)~elev1+elev2+elev3+elev4+maxNDVI1+seNDVI
                  +Matern(1|longitude+latitude), method="PQL", data=Loaloa,
                  family=binomial(), fixed=list(nu=0.5,rho=2.255197,lambda=1.075))   

 ## longer computation requiring interpolation of 197 points 
 ## Also illustrating effect of 'return.' argument
 res <- filled.mapMM(lfit,add.map=TRUE,plot.axes=quote({axis(1);axis(2)}),
    decorations=quote(points(pred[,coordinates],pch=15,cex=0.3)),
    return.="raster", # so that 'res' is a list representing a raster. 
    plot.title=title(main="Inferred prevalence, North Cameroon",
                     xlab="longitude",ylab="latitude"))
}

}
\keyword{hplot}
