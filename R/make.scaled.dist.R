make.scaled.dist <-
function(distMatrix,uniqueGeo,rho) {
  if ( missing(distMatrix) ) { ## 
     if ( missing(uniqueGeo) ) {
       mess <- pastefrom("missing(distMatrix) && missing(uniqueGeo).",prefix="(!) From ")
       stop(mess)
     } else {
        if (length(rho)==1) {
           uniqueScal<-uniqueGeo * rho 
        } else if (ncol(uniqueGeo)==length(rho)) {
           uniqueScal<-t(t(uniqueGeo) * rho) ## valid for vectorial rho...
        } else {
          mess <- pastefrom("invalid length(rho).",prefix="(!) From ")
          print(mess)
          mess  <- paste("Length should be either 1 or",ncol(uniqueGeo))
          stop(mess)
        }
     }
     ## scaled.dist <- as.matrix(dist(uniqueScal))
     scaled.dist <- dist(uniqueScal) ## one can do a lot with a dist object !
  } else if (length(rho)==1) {
     scaled.dist <- rho * distMatrix
  } else {
    mess <- pastefrom("input 'distMatrix' but length(rho)!=1.",prefix="(!) From ")
    stop(mess)
  }
  return(scaled.dist)
}
