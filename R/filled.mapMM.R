filled.mapMM <-
function(fitobject,coordinates,xrange=NULL,yrange=NULL,
                           margin=1/20,map.formula,phi=1e-05,gridSteps=41,
                           add.points=quote(points(pred[,coordinates],cex=1,lwd=2)),
                           add.map=FALSE,
#                           mask = 1,
                           axes = TRUE, plot.axes,map.asp=NULL,...) {
  if (missing(coordinates)) coordinates <- colnames(attr(fitobject,"info.uniqueGeo"))
  #coordinates <- c("longitude","latitude")  
  if (length(coordinates)!=2L) {
    stop(paste("'map' plots only 2D maps, while coordinates are of length ",length(coordinates),sep=""))
  }
  pred <- predict(fitobject)
  if (missing(map.formula)) map.formula <-  as.formula(paste("fitted ~ 1 + ",paste(findSpatial(fitobject$predictor))))
  smooobject <- corrHLfit(map.formula,data=pred,ranFix=list(phi=phi)) ## interpolates the predicted values...
  smoo <- predict(smooobject)
  x <- smoo[,coordinates[1]]
  y <- smoo[,coordinates[2]]
  if (is.null(xrange)) {
    xrange <- range(x)
    margex <- (xrange[2]-xrange[1])*margin
    xrange  <- xrange+margex*c(-1,1)
  }
  if (is.null(yrange)) {
    yrange <- range(y)
    margey <- (yrange[2]-yrange[1])*margin
    yrange  <- yrange+margey*c(-1,1)
  }  
  if (missing(plot.axes)) {
    if (axes) {
      title(main = "", xlab = "", ylab = "") ## ensures that nothing is draws and thus provides control to other parts of code
      Axis(x, side = 1)
      Axis(y, side = 2)
    }
  }
  
  xGrid <- seq(xrange[1],xrange[2],length.out=gridSteps)
  yGrid <- seq(yrange[1],yrange[2],length.out=gridSteps)
  newX <- expand.grid(xGrid,yGrid)
  colnames(newX) <- coordinates
  #  predictor <- smooobject$predictor
  #  yname <- paste(attr(predictor,"oriFormula"))[[2]]
  #  newX[[yname]] <- 0 
  gridpred <- predict(smooobject,newX = newX) # uses predictor of interpolation of predictor from orginal fit...
#  gridpred <- gridpred*blackcapMask #mask ## le test c'est blackcapMask
  if (is.logical(add.map)) {
    if(add.map) {
      if (require(maps)) {
        add.map <- quote(maps::map(,xlim=xrange,ylim=yrange,add=TRUE))
      } else {
        message("Package 'maps' not available, 'add.map' is ignored.")
        add.map <- quote(NULL)
      }  
    } else add.map <- quote(NULL) 
  } else add.map <- quote(add.map)
  spaMM.filled.contour(x=xGrid,y=yGrid,z=matrix(gridpred[,"fitted"],ncol=gridSteps),
                       plot.axes={plot.axes;
                                  eval(add.points);
                                  eval(add.map)
                       },map.asp=map.asp,...)  
}
