.get_qr <- function(mat,provide=TRUE) {
  envir <- attr(mat,"envir") ## environment
  envir$callcount <- envir$callcount+1L
  if (is.null(envir$qr)) { 
    if (provide) envir$qr <- .QRwrap_but_diag(mat)
  } 
  return(envir$qr)
}




