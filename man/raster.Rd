\name{raster}
\alias{raster}

\title{Raster plots of predicted responses}

\description{This documentation provides an example of producing a map out of spatial predictions, using rasters and geographical projections. No attempt has been made to provide a function automating the whole process, but the example can be used as a template.}
\note{This example was kindly provided by Alexandre Courtiol, Leibniz Instutite of Zoo and Wildlife Research, Berlin.}
\examples{
if (spaMM.getOption("example_maxtime")>10) {          
## using the Bonne projection for the fit *and* the plots

if(require("rgdal", quietly = TRUE) && require("rasterVis", quietly = TRUE)) {
  data(blackcap)
  
  ## fit:
  projection <- "+proj=bonne +lat_1=30n +lon_0=10e +units=km"
  bonne <- as.data.frame(project(as.matrix(blackcap[,c("longitude","latitude")]), projection))
  colnames(bonne) <- c("bonnex","bonney")
  bc <- cbind(blackcap,bonne)
  bfit <- corrHLfit(migStatus ~ 1 + Matern(1|bonnex+bonney), data=bc,
                    HLmethod="ML")
  
  ## plots:
  
  xrange <- range(bc$bonnex)
  yrange <- range(bc$bonney)
  gridSteps <- 101
  border <- 20
  xGrid <- seq(xrange[1]-border, xrange[2]+border, length.out = gridSteps)
  yGrid <- seq(yrange[1]-border, yrange[2]+border, length.out = gridSteps)
  newdata <- expand.grid(bonnex=xGrid, bonney=yGrid)
  gridpred <- predict(bfit, newdata = newdata)
  
  palette <- spaMM.colors()
  
  CreateRaster <- function(
    long,
    lat,
    values,
    proj="+proj=longlat +datum=WGS84",
    save.spatial.files=FALSE,
    filename="data_raster",
    overwrite.spatial.files=TRUE
  ) {
    # Args:
    #   long: a vector of the longitudes of the raster cells
    #   lat: a vector of the latitudes of the raster cells
    #   values: a vector of the values of the raster cells
    #   proj: the projection system for the raster
    #   save.spatial.files: logical indicating if 
    #          an hard copy of the raster should be saved (as ascii)
    #   filename: name of the file for the hard copy
    #   overwrite.spatial.files: logical indicating if 
    #          an existing hard copy should be overwritten or not
    #
    # Returns:
    #   The raster.
    #
    data <- data.frame(long=long, lat=lat, values=values)  # a dataframe 
    #                 with longitudes, lattitudes, and values is being created
    coordinates(data) <- ~long+lat  # coordinates are being set for the raster
    proj4string(data) <- CRS(proj)  # projection is being set for the raster
    gridded(data) <- TRUE  # a gridded structure is being set for the raster
    data.raster <- raster(data)  # the raster is being created
    if(save.spatial.files) writeRaster(
      data.raster,
      filename=paste(filename, ".asc", sep=""),
      overwrite=overwrite.spatial.files
    )  # if save=TRUE the raster is exported as an ascii file
    return(data.raster) # the raster is being returned
  }
  
  raster.predict <- CreateRaster(
    long=newdata$bonnex,
    lat=newdata$bonney,
    values=gridpred,
    proj=projection)
  
  data(worldcountries)
  data(oceanmask)
  # All shape files can be found here: http://www.naturalearthdata.com/downloads/
  # and loaded into R by e.g.
  # oceanmask <- readOGR("ne_110m_ocean.shp", layer="ne_110m_ocean")
  
  ## crop the ocean mask before projecting it, otherwise artefacts will appear 
  # (1) project world
  worldcountries <- spTransform(worldcountries, projection)
  # (2) define extent of region to be cropped in projected world
  extent.crop <- extent(projectExtent(raster.predict, oceanmask@proj4string))+20
  # (3) crop oceanmask according to this extent
  oceanmask <- crop(oceanmask, extent.crop)
  oceanmask <- spTransform(oceanmask, projection)
  plot(oceanmask)
  
  
  p <- levelplot( ## filled contour with legend bar
    raster.predict, ## levelplot from rasterVis handles any raster 
    maxpixels=4e6,
    margin=FALSE,
    cuts=length(palette)-1,
    col.regions=palette
  )
  
  
  country.layer <- layer(
    sp.polygons(worldcountries, fill=fill),
    data=list(sp.polygons=sp.polygons, worldcountries=worldcountries, 
              fill="transparent") 
  )
  
  p+country.layer
  
  CreateSpatialPoints <-
    function(
      long,
      lat,
      values=-9999,
      proj="+proj=longlat +datum=WGS84",
      save.spatial.files=FALSE,
      filename="my_points",
      overwrite.spatial.files=TRUE
    ) {
      data.sp <- data.frame(long=long, lat=lat, values=values)
      coordinates(data.sp) <- ~long+lat
      proj4string(data.sp) <- CRS(proj)
      if(save.spatial.files)   writeOGR(
        data.sp,
        dsn="./",
        layer=filename,
        morphToESRI=TRUE,
        driver="ESRI Shapefile",
        overwrite_layer=overwrite.spatial.files
      )
      return(data.sp)
    }
  
  my.points <- CreateSpatialPoints( ## userdefined
    long=bc$bonnex,
    lat=bc$bonney,
    proj=projection
  )
  
  points.layer <- layer(
    sp.points(my.points, col=col, cex=cex, pch=pch, lwd=lwd),
    data=list(sp.points=sp.points, my.points=my.points, 
              col="white", cex=2, pch=15, lwd=2)
  )
  
  p+points.layer
  
  ocean.layer <- layer(
    sp.polygons(oceanmask, fill=fill, col=col),
    data=list(sp.polygons=sp.polygons, oceanmask=oceanmask, 
              fill="black", col="transparent")
  )
  p+ocean.layer+country.layer+points.layer
  
}

}
}
