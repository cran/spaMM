is.square.diagonal <-
function( x, tol=1e-8 ) {
  if ( ! is.matrix( x ) ) stop( "argument x is not a numeric matrix" )
  if ( ncol(x) != nrow(x) ) return( FALSE )
  if ( ! is.numeric( x ) ) stop( "argument x is not a numeric matrix" )
  y <- x
  diag( y ) <- rep( 0, nrow( y ) )
  return( max(abs( y )) < tol )
}
