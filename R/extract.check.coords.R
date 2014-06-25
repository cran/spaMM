extract.check.coords <-
function(spatial.model,datanames) {
  if ( ! is.null(spatial.model)) {
    bars <- spatial.model[[2]] 
    coordinates <- DEPARSE(bars[[3]]) ## "x + y"
    coordinates <-  strsplit(coordinates," ")[[1]]
    coordinates <- setdiff(coordinates,c("+","%in%",":","/"))
  } else {
    stop("Call to 'HLCor' without a spatial term in the formula is suspect.")
    ## very old code handling old syntax with (1|pos) and default values of the coordinates argument
    coordinates <- c("x","y") ## back compat
  }
  coordcols <- which(datanames %in% coordinates)
  if ( length(coordcols) != length(coordinates) ) {
    stop("variables 'coordinates' not all found in the 'data'")
  }
  return(coordinates) ## should be ordered as bars[[3]] (important for predict)
}
