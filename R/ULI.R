ULI <-
function(...) {  ## Unique Location Index; '...' are simply names of variables in the data
 redondGeo<-cbind(...) ## always a matrix
 uniqueGeo<-unique(redondGeo) ## seems to preserve order
 ##FR->FR and I should in principle use a comparison of character (see ?unique)
 #### a 1 * nrow(redondGeo) matrix which nth element gives the index of the unique location in uniqueGeo :
 # designRU <- apply(redondGeo,1,function(v) {which(apply(v==t(uniqueGeo),2,all))}) ## this is slow; alternative using proxy:
 bla <- proxy::dist(uniqueGeo,redondGeo,method="Manhattan")
 designRU <- apply(bla,2,function(v){which(v==0L)})
 as.vector(designRU) ## has no names  
}
