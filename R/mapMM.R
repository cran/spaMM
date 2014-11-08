# il fallait fields pour migraine.colors: cf remplacement splint dans spaMM.colors

`spaMM.colors` <- function (n = 64, redshift = 1) {
  orig <- c("#00008F", "#00009F", "#0000AF", "#0000BF", "#0000CF", 
            "#0000DF", "#0000EF", "#0000FF", "#0010FF", "#0020FF", 
            "#0030FF", "#0040FF", "#0050FF", "#0060FF", "#0070FF", 
            "#0080FF", "#008FFF", "#009FFF", "#00AFFF", "#00BFFF", 
            "#00CFFF", "#00DFFF", "#00EFFF", "#00FFFF", "#10FFEF", 
            "#20FFDF", "#30FFCF", "#40FFBF", "#50FFAF", "#60FF9F", 
            "#70FF8F", "#80FF80", "#8FFF70", "#9FFF60", "#AFFF50", 
            "#BFFF40", "#CFFF30", "#DFFF20", "#EFFF10", "#FFFF00", 
            "#FFEF00", "#FFDF00", "#FFCF00", "#FFBF00", "#FFAF00", 
            "#FF9F00", "#FF8F00", "#FF8000", "#FF7000", "#FF6000", 
            "#FF5000", "#FF4000", "#FF3000", "#FF2000", "#FF1000", 
            "#FF0000", "#EF0000", "#DF0000", "#CF0000", "#BF0000", 
            "#AF0000", "#9F0000", "#8F0000", "#800000")
  orig[1:20] <- topo.colors(64)[1:20]
  if (n == 64 && redshift == 1) 
    return(orig)
  rgb.tim <- t(col2rgb(orig))
  temp <- matrix(NA, ncol = 3, nrow = n)
  x <- (seq(0, 1, , 64)^(redshift))
  xg <- seq(0, 1, , n)
  for (k in 1:3) {
    # hold <- splint(x, rgb.tim[, k], xg) ## using fields
    ## can be replaced by:
    colorpts <- data.frame(x=x,y=rgb.tim[,k])
    blob <- corrHLfit(y~ Matern(1|x),data=colorpts,ranFix=list(phi=1e-07,nu=4,rho=10)) ## use rho large...
    hold <- predict(blob,newX=data.frame(x=xg))[,"fitted"]
    ##
    hold[hold < 0] <- 0
    hold[hold > 255] <- 255
    temp[, k] <- round(hold)
  }
  rgb(temp[, 1], temp[, 2], temp[, 3], maxColorValue = 255)
}





`spaMM.filled.contour` <- function (x = seq(0, 1, length.out = nrow(z)), y = seq(0, 1, 
                                                       length.out = ncol(z)), z, xlim = range(x, finite = TRUE), 
          ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
          levels = pretty(zlim, nlevels), nlevels = 20, color.palette = spaMM.colors, 
          col = color.palette(length(levels) - 1), plot.title, plot.axes, 
          key.title, key.axes, map.asp = NULL, xaxs = "i", yaxs = "i", las = 1, 
          axes = TRUE, frame.plot = axes, ...) 
{
  if (missing(z)) {
    if (!missing(x)) {
      if (is.list(x)) {
        z <- x$z
        y <- x$y
        x <- x$x
      }
      else {
        z <- x
        x <- seq.int(0, 1, length.out = nrow(z))
      }
    }
    else stop("no 'z' matrix specified")
  }
  else if (is.list(x)) {
    y <- x$y
    x <- x$x
  }
  if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
    stop("increasing 'x' and 'y' values expected")
  mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
  on.exit(par(par.orig))

  xspan <- (xlim[2]-xlim[1])
  yspan <- (ylim[2]-ylim[1])
  if (is.null(map.asp)) map.asp <- yspan/xspan

  mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
  on.exit(par(par.orig))
  wscale <- (3 + mar.orig[2L]) * par("csi") * 2.54
  wmap <- par("din")[1]*2.54 - wscale
  Wmargin <- (par("din")[1]-par("pin")[1])*2.54
  wplotmap <- wmap - Wmargin  ## likely width of plot area
  Hmargin <- (par("din")[2]-par("pin")[2])*2.54
  hmap <- wplotmap*map.asp + Hmargin
  if (hmap>(par("din")[2]*2.54)) {
    hmap <- (par("din")[2]*2.54)
    reduction <- (hmap-Hmargin)/wplotmap
    wmap <- (hmap - Hmargin)/map.asp
  }
  layout(matrix(c(2, 1), ncol = 2L), widths = c(lcm(wmap),lcm(wscale)),heights=c(lcm(hmap)),respect=TRUE)
  
  par(las = las)
  mar <- mar.orig
  mar[4L] <- mar[2L]
  mar[2L] <- 1
  par(mar = mar)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i", 
              yaxs = "i")
  rect(0, levels[-length(levels)], 1, levels[-1L], col = col)
  if (missing(key.axes)) {
    if (axes) 
      axis(4)
  }
  else key.axes
  box()
  if (!missing(key.title)) 
    key.title
  mar <- mar.orig
  mar[4L] <- 1
  par(mar = mar)
  plot.new()
  plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs)
  .filled.contour(x, y, z, levels, col)
  if (missing(plot.axes)) {
    if (axes) {
      title(main = "", xlab = "", ylab = "")
      Axis(x, side = 1)
      Axis(y, side = 2)
    }
  }
  else plot.axes
  if (frame.plot) 
    box()
  if (missing(plot.title)) 
    title(...) ## not that this function does not know 'coordinates'
  else plot.title
  invisible()
} ## spaMM.filled.contour

mapMM <- function (fitobject,coordinates=NULL,xrange=NULL,yrange=NULL,margin=1/20,add.map= FALSE,
                   nlevels = 20, color.palette = spaMM.colors, 
                   map.asp=NULL,
                   col = color.palette(length(levels) - 1), 
                   plot.title, 
                   plot.axes, decorations,
                   key.title, key.axes, xaxs = "i", yaxs = "i", las = 1, 
                   axes = TRUE, frame.plot = axes,...) 
{
  
  if (missing(coordinates)) coordinates <- colnames(attr(fitobject,"info.uniqueGeo"))
  #coordinates <- c("longitude","latitude")  
  if (length(coordinates)!=2L) {
    stop(paste("'map' plots only 2D maps, while coordinates are of length ",length(coordinates),sep=""))
  }
  pred <- predict(fitobject)
  x <- pred[,coordinates[1]]
  y <- pred[,coordinates[2]]
  Zvalues <- pred[,"fitted"]
  zscaled<- 1+floor(nlevels*(0.000001+0.999998*(Zvalues-min(Zvalues))/(max(Zvalues)-min(Zvalues)))) ## makes sure its floor( ]1,nlevels+1[ ) 
  ZColor<- spaMM.colors(n=nlevels)  
  levels <- pretty(range(Zvalues), nlevels)
  if (is.null(xrange)) {
    xrange <- range(x)
  }
  xspan <- (xrange[2]-xrange[1])
  margex <- xspan * margin
  xrange  <- xrange+margex*c(-1,1)
  if (is.null(yrange)) {
    yrange <- range(y)
  }  
  yspan <- (yrange[2]-yrange[1])
  margey <- yspan * margin
  yrange  <- yrange+margey*c(-1,1)
  
  if (is.null(map.asp)) map.asp <- yspan/xspan
  
  mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
  on.exit(par(par.orig))
  wscale <- (3 + mar.orig[2L]) * par("csi") * 2.54
  wmap <- par("din")[1]*2.54 - wscale
  Wmargin <- (par("din")[1]-par("pin")[1])*2.54
  wplotmap <- wmap - Wmargin  ## likely width of plot area
  Hmargin <- (par("din")[2]-par("pin")[2])*2.54
  hmap <- wplotmap*map.asp + Hmargin
  if (hmap>(par("din")[2]*2.54)) {
    hmap <- (par("din")[2]*2.54)
    reduction <- (hmap-Hmargin)/wplotmap
    wmap <- (hmap - Hmargin)/map.asp
  }
  layout(matrix(c(2, 1), ncol = 2L), widths = c(lcm(wmap),lcm(wscale)),heights=c(lcm(hmap)),respect=TRUE)
  ## respect=TRUE est crucial...
#  layout(matrix(c(2, 1), ncol = 2L), widths = c(wmap,wscale),heights=c(hmap),respect=TRUE)
  par(las = las)
  mar <- mar.orig
  mar[4L] <- mar[2L]
  mar[2L] <- 1
  par(mar = mar)
  ## SCALE
  plot.new()
  zrange <- range(Zvalues)
  plot.window(xlim = c(0, 1), ylim = zrange,xaxs = "i", 
              yaxs = "i")
  rect(0, levels[-length(levels)], 1, levels[-1L], col = col)
  
  if (missing(key.axes)) {
    if (axes) 
      axis(4)
  }
  else key.axes
  box()
  if (!missing(key.title)) 
    key.title
  mar <- mar.orig
  mar[4L] <- 1
  par(mar = mar)
  #  plot.new()
  ## MAP
  plot.window(xrange, yrange, "", xaxs = xaxs, yaxs = yaxs)
  topontop <- order(zscaled,decreasing=FALSE) 
  plot(x=x[topontop],y=y[topontop], ## so that highest value will be printed last
       xlab="",ylab="", # 06/2015: give control to plot.title
       axes=FALSE, ## to retain control in later call
       col=ZColor[zscaled[topontop]],lwd=2)
  if (is.logical(add.map)) {
    if(add.map) {
      if (require(maps)) {
        maps::map(,xlim=xrange,ylim=yrange,add=TRUE)  ## require + :: is the way for objects from packages in Suggests:
      } else message("Package 'maps' not available, 'add.map' is ignored.")
    } 
  } else eval(add.map) ## the user may have included a map() in it but it's his problem...
  if (missing(plot.axes)) {
    if (axes) {
      title(main = "", xlab = "", ylab = "")
      Axis(x, side = 1)
      Axis(y, side = 2)
    }
  }
  else plot.axes
  if ( ! missing(decorations)) decorations
  if (frame.plot) box()
  dotlist <- list(...) 
  if (is.null(dotlist$xlab)) dotlist$xlab <- coordinates[1]
  if (is.null(dotlist$ylab)) dotlist$ylab <- coordinates[2]
  if (missing(plot.title)) 
    do.call(title,dotlist)
  else plot.title
  invisible()
} ## end mapMM

`filled.mapMM` <- function(fitobject,coordinates,xrange=NULL,yrange=NULL,
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
