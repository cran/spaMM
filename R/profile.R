spaMM.options <- function(..., warn=TRUE) {
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
  if (anyNA(names(old))) { # has element(s) $<NA>=NULL 
    if (warn) {
      checknames <- argnames[which(is.na(names(old)))]
      checknames <- setdiff(checknames, c("sparse_precision")) # exception for valid names not in default spaMM.options()
      if (length(checknames)) warning(paste0("'",paste(checknames,collapse="', '")),
                                      "' not previously in spaMM.options. Check such name(s)?", immediate. = TRUE)
    }
    names(old) <- argnames 
  }
  .spaMM.data$options[argnames] <- temp
  invisible(old)
}

spaMM.getOption <- function (x) {spaMM.options(x, warn=FALSE)[[1]]}


## large rho is not a problem
## large nu is a problem the more so when rho is small (=> 'small scaled distance gaussian')
# lme did not manage to go beyond nu=17.48 in case Vn phiFix...

".onAttach" <- function (lib, pkg) {
  version <- utils::packageVersion("spaMM")
  packageStartupMessage("spaMM (Rousset & Ferdy, 2014, version ", version, 
                          ## not sure this will always work and makes sense only for devel version :
                          # ", packaged ", utils::packageDescription("spaMM")$Packaged,
                        ") is loaded.", 
    "\nType 'help(spaMM)' for a short introduction,",
    "\n'news(package='spaMM')' for news,",
    "\nand 'citation('spaMM')' for proper citation.")
  #unlockBinding(".SpaMM", asNamespace("spaMM")) ## required when a .SpaMM list was used instead of an envir
  
}

.onLoad <- function(libname, pkgname) {
  if ( ! proxy::pr_DB$entry_exists("Earth")) {
    pr_DB$set_entry(FUN = .Dist.earth.mat, names = c("Earth", "dist.earth"))
    pr_DB$modify_entry(
      names = "Earth",
      description = "Approximate great-circle distance in Km between points on Earth surface.",
      loop = FALSE,
      distance = TRUE
    )
  } else warning("'Earth' entry already present in proxy::pr_DB database.")
  if ( ! proxy::pr_DB$entry_exists("EarthChord")) {
    pr_DB$set_entry(FUN = .Dist.chord.mat, names = c("EarthChord", "dist.EarthChord"))
    pr_DB$modify_entry(
      names = "EarthChord",
      description = "Approximate chord distance in Km between points on Earth surface.",
      loop = FALSE,
      distance = TRUE
    )
  } else warning("'EarthChord' entry already present in proxy::pr_DB database.")
  success <- suppressMessages(do.call("require",list(package="memoise", quietly=TRUE))) # 'quietly' needed to suppress *warning* when memoise is not attached.
  .spaMM.data$options$need_memoise_warning <- ! success
  if (success) { # remarkably, no need for special wrapper. R CMD check handling of .onLoad must be special
    ..CMP_mu2lambda <<- memoise(f=..CMP_mu2lambda, omit_args="CMP_linkfun_objfn", 
                                cache = .do_call_wrap("cache_mem", arglist=list(max_size = 10 * 1024^2), pack="cachem"))
    .Rcpp_COMP_Z <<- memoise(f=.Rcpp_COMP_Z, cache = .do_call_wrap("cache_mem", arglist=list(max_size = 10 * 1024^2), pack="cachem"))
    #  str(environment(spaMM:::.Rcpp_COMP_Z)$"_cache"$keys()) to get info on the cache...
    ..trDiagonal <<- memoise(f=..trDiagonal, cache = .do_call_wrap("cache_mem", arglist=list(max_size = 1024^2), pack="cachem"))
    .get_phantom_map <<- memoise(f=.get_phantom_map, cache = .do_call_wrap("cache_mem", arglist=list(max_size = 1024^2), pack="cachem"))
  } 
}

".onUnload" <- function (libpath) {
  pr_DB$delete_entry("Earth")
  pr_DB$delete_entry("EarthChord")
  library.dynam.unload("spaMM", libpath)
} ## testable by calling unloadNamespace("spaMM")
#  pkgpath <- system.file(package="OKsmooth") # https://github.com/hadley/devtools/issues/119

# unloadNampespace() calls .onUnload only after after checking dependencies, so the following would be useless in .onUnload()
.unloads4spaMM <- function() {
  unloadNamespace("probitgem")
  unloadNamespace("IsoriX")
  unloadNamespace("Infusion")
  unloadNamespace("blackbox")
}

.Dist.earth.mat <- function (x, y=NULL, radius=6371.009) { # x and y are both matrices. In each, first col is longitude, second is latitude
  ## Earth radius used for approximation = 6371.009 = 1/3*(2*6378.137+6356.752)  [details on https://en.wikipedia.org/wiki/Great-circle_distance]
  ## This function computes orthodromic distances in Km between locations.
  rad_deg <- pi/180 # input should be in degrees, and this converts to radians.
  x <- x*rad_deg
  if(is.null(y)) { ## distances within matrix
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
  pp <- radius * acos(pp)
  if (is.null(y)) pp <- as.dist(pp)  ## spaMM wants an half matrix in this case, not a full one
  return(pp)
}

.Dist.chord.mat <- function(x,y=NULL, radius=6371.009) { # part of EarthChord implementation
  pp <- .Dist.earth.mat(x,y,radius=1)
  pp <- radius * 2*sin(pp/2)
  if (is.null(y)) pp <- as.dist(pp)  ## spaMM wants an half matrix in this case, not a full one
  return(pp)
}
