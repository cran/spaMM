if (FALSE) { ## version de developpement a conserver
  .create_get_qr <- function(tag) {
    ## creates the function in a local envir different from the closure of HLfit
    # so that this envir can be manipulated easily
    qr <- NULL## so that no build/check note on the <<- 
    callcount <- NULL
    locfn <- function(mat,provide=TRUE) {
      debugtag <- "wAugX"
      callcount <<- callcount+1L
      if (is.null(qr)) { 
        ## This call gets args (except res) from the envir of locfn def'd below:
        if (provide) {
          qr <<- QRwrap(mat)
          if (tag==debugtag) {cat(paste(" provide ",tag))}
        } else if (tag==debugtag) {cat(paste(" NULL",tag)); if (callcount>1L) browser()} ## can stop at end of HLfit... 
      } else if (tag==debugtag) cat(paste(" OK",tag))
      return(qr)
    } ## this function changes environment(<res object>$get_info_crits)$info_crits
    environment(locfn) <- list2env(list(qr=NULL,tag=tag,callcount=0L),
                                   parent=environment(QRwrap))
    return(locfn)
  }
  #mat <- structure(matrix(runif(12),ncol=3),get_qr=create_get_qr())
  #essai <- function(mat) {attr(mat,"get_qr")(mat)}
  #essai(mat) ##changes mat 'globally'
}  else {
  .get_qr <- function(mat,provide=TRUE) {
    envir <- attr(mat,"envir") ## environment
    envir$callcount <- envir$callcount+1L
    if (is.null(envir$qr)) { 
      if (provide) envir$qr <- QRwrap(mat)
    } 
    return(envir$qr)
  }
}



