Variogram.corMatern <-
function(object, distance = NULL, sig2 = 1, length.out = 50, FUN, ...)
{
  if (is.null(distance)) {
    rangeDist <- range(unlist(getCovariate(object)))
    distance <- seq(rangeDist[1], rangeDist[2], length = length.out)
  }
  params <- coef(object, unconstrained = FALSE)
  if (length(params) == 2) {            # no nugget effect
    range  <- params[1]
    nu   <- params[2]
    nugg <- 0
  } else {                              # nugget effect
    range  <- params[1]
    nu   <- params[2]
    nugg <- params[3]
  }
  cat(range,nu,nugg,"\n")
  f <- function(distance) {
    z <- rep(0,length(distance))
    for(i in 1:length(distance)) {
      z[i] <- .C("matern_cor",
         c(range,nu,nugg),
         as.double(distance[i]),
         as.integer(attr(object, "nugget")),
         as.integer(attr(object, "nuScaled")),
         double(1))[[5]]
    }
    return(z)
  }
  ## static void matern(double *par, double dist, longint *nug, longint *nusc,
  ##          double *cor)

  zz <- f(distance)
  print(zz)
  val <- data.frame(variog = sig2 * (nugg + (1 - nugg) * zz),
                    dist = distance)
  class(val) <- c("Variogram", "data.frame")
  val
}
