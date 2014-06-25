post.process.ZALlist <-
function(ZAL,predictor) {
  nrand <- length(ZAL)
  if (nrand>1) {
    ZAL <- do.call(cbind,ZAL)
    # FR->FR ici ça serait bien de pouvoir identifier un diagonal matrix...
  } else if (nrand==1) {
    ZAL <- ZAL[[1]] 
    if ((! is.null(attr(predictor,"%in%"))) && attr(predictor,"%in%") && ncol(ZAL)==nrow(ZAL)) { 
      ## test of the attribute is a heuristic way of detecting when using the block structure will lead to faster analysis
      partition <- findblocks(ZAL) 
      if ( length(partition)>1 ) {
        partition <- cumsum(c(0,partition))
        attr(ZAL,"partition") <- partition
      }
    }
  }
  return(ZAL)
}
