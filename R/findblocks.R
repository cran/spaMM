findblocks <-
function(mat,symm=FALSE) {
  if (!symm) {
    mat <- abs(mat)
    mat <- mat+t(mat)
  }
  nc <- ncol(mat)
  if (nc>1) {
    for (it in (2:nc)) {
      chk <- max(mat[1:(it-1),it:nc,drop=FALSE])
      if(chk==0) break
    }  
    if (it==nc) {
      return(nc)
    } else {
      return(c(it-1,findblocks(mat[it:nc,it:nc,drop=FALSE],symm=TRUE)))
    }
  } else {
    return(1)
  }
}
