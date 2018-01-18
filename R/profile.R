spaMM.options <- function(...) {
  if (nargs() == 0) return(.spaMM.data$options)
  temp <- list(...)
  if (length(temp) == 1 && is.null(names(temp))) {
    arg <- temp[[1]]
    switch(mode(arg),
           list = temp <- arg,
           character = return(.spaMM.data$options[arg]),  ## return here for eg ... = "NUMAX"
           stop("invalid argument: ", sQuote(arg)))
  }
  if (length(temp) == 0) return(.spaMM.data$options)
  argnames <- names(temp)
  if (is.null(argnames)) stop("options must be given by name")
  old <- .spaMM.data$options[argnames]
  names(old) <- argnames ## bc names are not valid for previously absent elements
  .spaMM.data$options[argnames] <- temp
  invisible(old)
}

spaMM.getOption <- function (x) {spaMM.options(x)[[1]]}


## large rho is not a problem
## large nu is a problem the more so when rho is small (=> 'small scaled distance gaussian')
# lme did not manage to go beyond nu=17.48 in case Vn phiFix...

".onAttach" <- function (lib, pkg) {
  version <- utils::packageVersion("spaMM")
  packageStartupMessage("spaMM (version ", version, 
                          ## not sure this will always work and makes sense only for devel version :
                          # ", packaged ", utils::packageDescription("spaMM")$Packaged,
                        ") is loaded.", 
    "\nType 'help(spaMM)' for a short introduction,\nand news(package='spaMM') for news.")
  #unlockBinding(".SpaMM", asNamespace("spaMM")) ## required when a .SpaMM list was used instead of an envir
  
}

.onLoad <- function(libname, pkgname) {
  ## Let us set the new proxy method
  if ( ! proxy::pr_DB$entry_exists("Earth")) {
    pr_DB$set_entry(FUN = .Dist.earth.mat, names = c("Earth", "dist.earth"))
    pr_DB$modify_entry(
      names = "Earth",
      description = "Approximate distance in Km between points on earth surface.",
      loop = FALSE,
      distance = TRUE
    )
  }
}

".onUnload" <- function (libpath) {
  pr_DB$delete_entry("Earth")
  library.dynam.unload("spaMM", libpath)
} ## testable by calling unloadNamespace("spaMM")
#  pkgpath <- system.file(package="OKsmooth") # https://github.com/hadley/devtools/issues/119

.Dist.earth.mat <- function (x, y=NULL) { # x and y are both matrices. In each, first col is longitude, second is latitude
  ## This function computes orthodromic distances in Km between locations.
  rad_deg <- pi/180
  x <- x*rad_deg
  if(is.null(y)) { ## distances within matrice
    coslat <- cos(x[, 2]) ## [,2] is latitude
    sinlat <- sin(x[, 2])
    coslon <- cos(x[, 1]) ## [,1] is longitude
    sinlon <- sin(x[, 1])
    pp <- cbind(coslat * coslon, coslat * sinlon, sinlat) %*% 
      t(cbind(coslat * coslon, coslat * sinlon, sinlat))
  } else { ## cross-matrices distances
    y <- y*rad_deg
    coslat1 <- cos(x[, 2])
    sinlat1 <- sin(x[, 2])
    coslon1 <- cos(x[, 1])
    sinlon1 <- sin(x[, 1])
    coslat2 <- cos(y[, 2])
    sinlat2 <- sin(y[, 2])
    coslon2 <- cos(y[, 1])
    sinlon2 <- sin(y[, 1])
    pp <- cbind(coslat1 * coslon1, coslat1 * sinlon1, sinlat1) %*% 
      t(cbind(coslat2 * coslon2, coslat2 * sinlon2, sinlat2))
  }
  pp <- pmin(pmax(pp,-1),1)
  ## Earth radius used for approximation = 6371.009 = 1/3*(2*6378.137+6356.752)  [details on https://en.wikipedia.org/wiki/Great-circle_distance]
  pp <- 6371.009 * acos(pp)
  if (is.null(y)) pp <- as.dist(pp)  ## spaMM wants an half matrix in this case, not a full one
  return(pp)
}
