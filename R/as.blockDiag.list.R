as.blockDiag.list <-
function(matlist,fill=0L) { ## FR->FR but no check that the matrix is square 
  bd <- matlist
  attr(bd,"rowNames") <- unlist(lapply(bd,rownames)) ## may be NULL
  attr(bd,"cumsum") <- cumsum(c(0,unlist(lapply(bd,ncol)))) 
  attr(bd,"fill") <- fill
  class(bd) <- c("blockDiag",class(bd))
  bd  
}
