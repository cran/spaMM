solveWrap.vector <-
function(A,b,...) {
  if (inherits(A,"diagonalMatrix")) return(b/diag(A))
  if (inherits(A,"Rcpp-QR")) { ## no pivoting, $R is upper triangular
    b <- t(t(b) %*% (A$Q))
    solved <- try(backsolve(A$R,b),silent=TRUE)
    return(solved) ## gives control to calling function 
  }
  ## all other cases
  safesolve.qr.vector(A,b,...)
}
