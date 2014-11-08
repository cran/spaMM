is.square.diagonal <- function( x, tol=1e-8 ) {
  if ( ! is.matrix( x ) ) stop( "argument x is not a numeric matrix" )
  if ( ncol(x) != nrow(x) ) return( FALSE )
  if ( ! is.numeric( x ) ) stop( "argument x is not a numeric matrix" )
  y <- x
  diag( y ) <- rep( 0, nrow( y ) )
  return( max(abs( y )) < tol )
}

is.identity <- function( x, tol=1e-8 ) {
  if (!is.square.diagonal( x ) ) return( FALSE )
  return(all(abs(range(diag(x))-1)<tol))
}

