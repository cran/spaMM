make.scaled.dist <-
function(distMatrix,uniqueGeo,rho,rho.mapping=seq_len(length(rho))) {
  if ( missing(distMatrix) ) { ## 
     if ( missing(uniqueGeo) ) {
       mess <- pastefrom("missing(distMatrix) && missing(uniqueGeo).",prefix="(!) From ")
       stop(mess)
     } else {
        if (length(rho)==1) {
           uniqueScal<-uniqueGeo * rho 
        } else if (ncol(uniqueGeo)==length(rho.mapping)) {
           uniqueScal<-t(t(uniqueGeo) * rho[rho.mapping]) ## valid for vectorial rho...
        } else {
          mess <- pastefrom("invalid length(rho[rho.mapping]).",prefix="(!) From ")
          print(mess)
          mess  <- paste("Length should be either 1 or",ncol(uniqueGeo))
          stop(mess)
        }
     }
     scaled.dist <- proxy::dist(uniqueScal) ## one can do a lot with a dist object !
     ## here scaled.dist is always a dist object: if uniqueScal was a single row, it is dist(0) 
  } else if (length(rho)==1) {
     scaled.dist <- rho * distMatrix
     ## here inconsistent behavior: scaled.dist is a dist object except if distMatrix was dist(0), scaled.dist is now numeric(0); we standardise it:
     if ( identical(scaled.dist, numeric(0))) scaled.dist <- dist(0)
  } else {
    mess <- pastefrom("input 'distMatrix' but length(rho)!=1.",prefix="(!) From ")
    stop(mess)
  }
  return(scaled.dist)
}
