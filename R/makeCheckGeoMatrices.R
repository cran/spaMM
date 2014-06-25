makeCheckGeoMatrices <-
function(data,distMatrix=NULL,uniqueGeo=NULL,coordinates) {
  isListData <- inherits(data,"list")
  if (is.null(distMatrix)) { 
    if ( is.null(uniqueGeo) ) { ## then construct it from the data ## this should be the routine case
      if (isListData) {
        uniqueGeo <- lapply(data,function(dd) {unique(dd[,coordinates,drop=FALSE])})
        nbUnique <- lapply(uniqueGeo,nrow) 
      } else {
        uniqueGeo <- unique(data[,coordinates,drop=FALSE])
        nbUnique <- nrow(uniqueGeo) 
      }
    } 
    ## (2): we need distMatrix *here* in all cases for the check
    if (isListData) {
      distMatrix <- lapply(uniqueGeo,proxy::dist)
    } else distMatrix <- proxy::dist(uniqueGeo)
  } else { ## there is a distMatrix, this is what will be used by HLCor
    if (isListData) {
      nbUnique <- lapply(seq_len(length(data)),function(dd) {checkDistMatrix(distMatrix,dd,coordinates)})
    } else nbUnique <- checkDistMatrix(distMatrix,data,coordinates)
    ## stops if problems, otherwise checkDistMatrix has no useful return value
  }
  return(list(nbUnique=nbUnique,uniqueGeo=uniqueGeo,distMatrix=distMatrix))
}
