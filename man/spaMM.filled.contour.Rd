\name{spaMM.filled.contour}
\alias{spaMM.filled.contour}

\title{Level (Contour) Plots with better aspect ratio control (for geographical maps, at least)}
\usage{
spaMM.filled.contour(x = seq(0, 1, length.out = nrow(z)),
               y = seq(0, 1, length.out = ncol(z)),
               z,
               xrange = range(x, finite = TRUE),
               yrange = range(y, finite = TRUE),
               zrange = range(z, finite = TRUE),
               margin=1/20,
               levels = pretty(zrange, nlevels), nlevels = 20,
               color.palette = spaMM.colors,
               col = color.palette(length(levels) - 1),
               plot.title, plot.axes, key.title=NULL, key.axes=NULL,
               map.asp = NULL, xaxs = "i", yaxs = "i", las = 1,
               axes = TRUE, frame.plot = axes, \dots)

}
\arguments{
 \item{x, y}{locations of grid lines at which the values in \code{z} are
   measured.  These must be in ascending order.  (The rest of this
   description does not apply to \code{.filled.contour}.)
   By default, equally spaced values from 0 to 1 are used.  If \code{x}
   is a \code{list}, its components \code{x$x} and \code{x$y} are used
   for \code{x} and \code{y}, respectively.  If the list has component
   \code{z} this is used for \code{z}.}
 \item{z}{a numeric matrix containing the values to be plotted..  Note that
   \code{x} can be used instead of \code{z} for convenience.}
 \item{xrange}{x range of the plot.}
 \item{yrange}{y range of the plot.}
 \item{zrange}{z range of the plot.}
\item{margin}{
   This controls how far (in relative terms) the plot extends beyond the x and y ranges 
   of the analyzed points, and is overriden by explicit \code{xrange} and \code{yrange} arguments.  
 }
 \item{levels}{a set of levels which are used to partition the range
    of \code{z}.  Must be \bold{strictly} increasing (and finite).  Areas
    with \code{z} values between consecutive levels are painted with the
    same color.}
 \item{nlevels}{if \code{levels} is not specified, the range of \code{z},
    values is divided into approximately this many levels.}
 \item{color.palette}{a color palette function to be used to assign
    colors in the plot.}
 \item{col}{an explicit set of colors to be used in the plot.
    This argument overrides any palette function specification.  There
    should be one less color than levels}
 \item{plot.title}{statements which add titles to the main plot.}
 \item{plot.axes}{statements which draw axes (and a \code{\link{box}})
   on the main plot.  This overrides the default axes.}
 \item{key.title}{statements which add titles for the plot key.}
 \item{key.axes}{statements which draw axes on the plot key.
    This overrides the default axis.}
\item{map.asp}{  
the y/x aspect ratio of the 2D plot area (not of the full figure including the scale). Default is (plotted y range)/(plotted x range) (i.e., scales for x are identical). 
}
 \item{xaxs}{the x axis style.  The default is to use internal
    labeling.}
 \item{yaxs}{the y axis style.  The default is to use internal
    labeling.}
 \item{las}{the style of labeling to be used.  The default is to
    use horizontal labeling.}
 \item{axes, frame.plot}{logicals indicating if axes and a box should be
   drawn, as in \code{\link{plot.default}}.}
 \item{\dots}{additional \link{graphical parameters}, currently only passed to
   \code{\link{title}()}.}
}
\description{
  This function is derived from \code{filled.contour} in the \code{graphics} package, and 
  this documentation is likewise heavily based on that of \code{filled.contour}.

  This function likewise produces a contour plot with the areas between the
  contours filled in solid color, and a 
  key showing how the colors map to z values is likewise shown to the right of
  the plot.
  The only difference is the way the aspect ratio is determined and can be controlled (using the \code{map.asp} parameter instead of \code{asp}),
  They thus easily provide nice-looking maps with meaningful latitude/longitude ratio (see Examples). However, this does not work well with rstudio.
}
\details{
  The values to be plotted can contain \code{NA}s.  Rectangles with two
  or more corner values are \code{NA} are omitted entirely: where there
  is a single \code{NA} value the triangle opposite the \code{NA} is
  omitted.

  Values to be plotted can be infinite: the effect is similar to that
  described for \code{NA} values.

}
\references{
  Cleveland, W. S. (1993)
  \emph{Visualizing Data}.
  Summit, New Jersey: Hobart.
}
\note{Builds heavily on filled.contour by Ross Ihaka and R-core.
  \code{spaMM.filled.contour} uses the \code{\link{layout}} function and so is
  restricted to a full page display.

  The output produced by \code{spaMM.filled.contour} is actually a combination
  of two plots; one is the filled contour and one is the legend.  Two
  separate coordinate systems are set up for these two plots, but they
  are only used internally -- once the function has returned these
  coordinate systems are lost.  If you want to annotate the main contour
  plot, for example to add points, you can specify graphics commands in
  the \code{plot.axes} argument.  See the Examples.
}
\seealso{
  \code{\link{contour}}, \code{\link{image}},
  \code{\link{palette}}; \code{\link[lattice:levelplot]{contourplot}}
  and \code{\link[lattice]{levelplot}} from package \code{lattice}.
}
\examples{
spaMM.filled.contour(volcano, color = spaMM.colors) # simple

## Comparing the layout with that of filled.contour:
#  (except that it does not always achieve the intended effect 
#  in RStudio Plots pane). 

x <- 10*1:nrow(volcano)
y <- 10*1:ncol(volcano)
spaMM.filled.contour(x, y, volcano, color = terrain.colors,
    plot.title = title(main = "The Topography of Maunga Whau",
    xlab = "Meters North", ylab = "Meters West"),
    plot.axes = { axis(1, seq(100, 800, by = 100))
                  axis(2, seq(100, 600, by = 100)) },
    key.title = title(main = "Height\n(meters)"),
    key.axes = axis(4, seq(90, 190, by = 10)))  # maybe also asp = 1
mtext(paste("spaMM.filled.contour(.) from", R.version.string),
      side = 1, line = 4, adj = 1, cex = .66)

## compare with      

filled.contour(x, y, volcano, color = terrain.colors,
    plot.title = title(main = "The Topography of Maunga Whau",
    xlab = "Meters North", ylab = "Meters West"),
    plot.axes = { axis(1, seq(100, 800, by = 100))
                  axis(2, seq(100, 600, by = 100)) },
    key.title = title(main = "Height\n(meters)"),
    key.axes = axis(4, seq(90, 190, by = 10)))  # maybe also asp = 1
mtext(paste("filled.contour(.) from", R.version.string),
      side = 1, line = 4, adj = 1, cex = .66)

}
