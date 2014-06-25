is.identity <-
function( x, tol=1e-8 ) {
  if (!is.square.diagonal( x ) ) return( FALSE )
  return(all(abs(range(diag(x))-1)<tol))
}
