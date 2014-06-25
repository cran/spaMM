solveWrap.matrix <-
function(A,B,...) {
  if (inherits(A,"diagonalMatrix")) return(B/diag(A)) ## works if A is the matrix, not its diagonal...
  if (inherits(A,"Rcpp-QR")) { ## no pivoting, $R is upper triangular
    solved <- try(backsolve(A$R,t(A$Q) %*%B),silent=TRUE)
    return(solved) ## gives control to calling function 
  }
  ## all other cases
  safesolve.qr.matrix(A,B,...)
}
