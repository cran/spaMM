make.scaled.dist <-
function(distm,uniqueGeo,rho) {
  if ( missing(distm) ) { ## 
     if ( missing(uniqueGeo) ) {
       mess <- pastefrom("missing(distm) && missing(uniqueGeo).",prefix="(!) From ")
       stop(mess)
     } else {
        if (length(rho)==1) {
           uniqueScal<-uniqueGeo * rho 
        } else if (ncol(uniqueGeo)==length(rho)) {
           uniqueScal<-t(t(uniqueGeo) * rho) ## valid for vectorial rho...
        } else {
          mess <- pastefrom("invalid length(rho).",prefix="(!) From ")
          stop(mess)
        }
     }
     ## scaled.dist <- as.matrix(dist(uniqueScal))
     scaled.dist <- dist(uniqueScal) ## one can do a lot with a dist object !
  } else if (length(rho)==1) {
     scaled.dist <- rho * distm
  } else {
    mess <- pastefrom("input 'distm' but length(rho)!=1.",prefix="(!) From ")
    stop(mess)
  }
  return(scaled.dist)
}
