ULI <-
function(...) {  ## Unique Location Index; '...' are simply names of variables in the data
 redondGeo<-cbind(...) ## always a matrix
 uniqueGeo<-unique(redondGeo) ## seems to preserve order
 #### a 1 * nrow(redondGeo) matrix which nth element gives the index of the unique location in uniqueGeo :
 designRU <- apply(redondGeo,1,function(v) {which(apply(v==t(uniqueGeo),2,all))}) ## has no names
 ##FR->FR and I should in principle use a comparison of character (see ?unique)
 as.vector(designRU)  
}
