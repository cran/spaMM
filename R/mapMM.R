mapMM <-
function (fitobject,coordinates=NULL,xrange=NULL,yrange=NULL,margin=1/20,add.map= FALSE,
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
}
