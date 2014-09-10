`%id*id%` <-
function(A,B) {
  if ( ( ! is.null(attr(A,"identityMatrix"))) && attr(A,"identityMatrix")) return(B)
  if ( ( ! is.null(attr(B,"identityMatrix"))) && attr(B,"identityMatrix")) return(A)
  return(A %*% B)
}
