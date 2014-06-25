checkDistMatrix <-
function(distMatrix,data,coordinates) {
  classDistm <- class(distMatrix)
  if ( ! classDistm %in% c("matrix","dist")) {
    message(paste("(!) 'distMatrix' argument appears to be a '",classDistm,"',",sep=""))
    stop("not a 'matrix' or 'dist'. Check the input. I exit.")
  }
  ## chol() fails on distances matrices with repeated locations (which are pos SD)... but chol() not used by default
  ## the following code assumes that distMatrix deals only with unique locations, and checks this
  ## HENCE ******* distMatrix must refer to unique values of a grouping variable *********
  usernames <- rownames(distMatrix)
  checknames <- all(sapply(usernames,function(v) {v %in% rownames(data)})) ## 
  if (!checknames) {
    warning("The rownames of 'distMatrix' are not rownames of the 'data'. Further checking of 'distMatrix' is not possible.")
    nbUnique <- NA
  } else {
    uniqueGeo <- unique(data[usernames,coordinates,drop=FALSE]) ## check that this corresponds to unique locations
    nbUnique <- nrow(uniqueGeo)
    if (nbUnique != nrow(distMatrix)) {
      stop("The dimension of 'distMatrix' does not match the number of levels of the grouping variable")
    } else { ## check order
      redondGeo <- data[,coordinates,drop=F]
      designRU <- apply(redondGeo,1,function(v) {which(apply(v==t(uniqueGeo),2,all))}) ## has no names
      ## eg 1 1 2 2 3 2 3 4 is valid for 8 obs, 4 unique locations
      designRU <- unique(as.vector(designRU)) ## should then be 1 2 3 4
      ## but if distMatrix in reverse order, the first row of redondGeo would match the 4th of uniqueGeo and then the following test is FALSE:
      if ( ! all (designRU==seq_len(length(designRU))) ) {
        stop("The rows of 'distMatrix' are not ordered as rows of the 'data'.")
      }
    } 
  }
  nbUnique ## if stop() did not occur
}
