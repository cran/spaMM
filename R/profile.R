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


if (FALSE) {
  # currently triggers 
  # "error: there is no package called 'rlang'" when testing on CRAN
  # if installed package can be loaded from temporary location.
  .cli_packageStartupMessage <- function() {
    mess <- paste0("spaMM (Rousset & Ferdy, 2014, version ", version, 
                   ## not sure this will always work and makes sense only for devel version :
                   # ", packaged ", utils::packageDescription("spaMM")$Packaged,
                   ") is loaded.", 
                   "\nSee {.topic [spaMM](spaMM::spaMM)} for a short introduction,",
                   "\n'news(package='spaMM')' for news,",
                   "\nand 'citation('spaMM')' for proper citation.",
                   "\nFurther infos, slides, etc. at https://gitlab.mbb.univ-montp2.fr/francois/spamm-ref.\n")
    # https://github.com/r-lib/cli/issues/589:
    cli::cli_inform(mess, class = "packageStartupMessage")
    if (.spaMM.data$options$dec2spp) 
      packageStartupMessage(cli::style_bold(cli::style_underline("This development version uses 'spprec' method\n")),
                            "in many cases where 'decorr' one has been previously used.")
  }
}

# additional (wrt .onLoad) operations when the package is visible to the user (:: not required to call a function)
".onAttach" <- function (lib, pkg) {
  version <- utils::packageVersion("spaMM")
  # if (FALSE) {
  #   .cli_packageStartupMessage() 
  # } else {
    mess <- paste0("spaMM (Rousset & Ferdy, 2014, version ", version, 
                   ## not sure this will always work and makes sense only for devel version :
                   # ", packaged ", utils::packageDescription("spaMM")$Packaged,
                   ") is loaded.", 
                   "\nSee 'help('spaMM')' for a short introduction,",
                   "\n'news(package='spaMM')' for news,",
                   "\nand 'citation('spaMM')' for proper citation.",
                   "\nFurther infos, slides, etc. at https://gitlab.mbb.univ-montp2.fr/francois/spamm-ref.\n")
    packageStartupMessage(mess)
  # }
}

# Whatever's needed for operation of the namespace (allows ::: or ::)
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
  backports::import(pkgname, "...names") # to ensure back compat as long as spaMM supports R < 4.1
  .setNbThreads(thr=1L) # at C++ level for Eigen; only initialization, can be modified by control.HLfit$nbTHreads
  .spaMM.data$options$Matrix_old <- (packageVersion("Matrix")<"1.4-2") # 1st Matrix public version 1.5-0 (2022-09-09 r3636)
  # so if I wait 4 years to tidy Matrix_old: september 2026...
  # Version 1.7 of Matrix requires R 4.4.0 so if I require it I change my R requirement
  # earlier Matrix version required R 3.5.0, but this was a bit buggy before Matrix version 1.6-2 (2023-11-05 r4503)
  # => require version 1.6-2 and don't wait 4 years ?
  .spaMM.data$options$HLnames <- unique(c(names(formals(HLCor)),names(formals(HLfit)), 
                                   "ADFun", # so that this private arg, in the dots, causes no warning and is passed to .preprocess() 
                                   names(formals(mat_sqrt)),names(formals(make_scaled_dist))))
}

".onUnload" <- function (libpath) {
  pr_DB$delete_entry("Earth")
  pr_DB$delete_entry("EarthChord")
  library.dynam.unload("spaMM", libpath)
} ## testable by calling unloadNamespace("spaMM")

# In contexts where libpath is not available, a syntax is 
#pd.file <- attr(packageDescription("spaMM"), "file")
#library.dynam.unload("spaMM", libpath = sub("/Meta.*", '', pd.file))

# unloadNampespace() calls .onUnload only after after checking dependencies, so the following would be useless in .onUnload()
.unloads4spaMM <- function() {
  unloadNamespace("probitgem")
  unloadNamespace("IsoriX")
  unloadNamespace("gspace2infr")
  unloadNamespace("Infusion")
  unloadNamespace("blackbox")
  unloadNamespace("spaMM") 
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
