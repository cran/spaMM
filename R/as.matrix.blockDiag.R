as.matrix.blockDiag <-
function(x,fill=attr(x,"fill"),...) {
  rowNames <- attr(x,"rowNames")
  Ncol <- length(rowNames)
  value <- array(fill, c(Ncol, Ncol), list(rowNames,rowNames)) ## full matrix of 0's with rows and col names
  for (i in seq_along(x)) {
    aux <- as.matrix(x[[i]],...) ## recursif...
    nam <- rownames(aux)
    value[nam,nam] <- as.vector(aux) ## copie du block dans le block indexé par ses col et row names
  }
  rownames(value) <- colnames(value) <- rowNames
  value
}
