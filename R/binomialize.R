binomialize <-
function(data,responses,sortedTypes=NULL,binResponse=c("npos","nneg"),depth=Inf,input="types") {
  if (input=="types") {
    counts <- table(unlist(data[,responses]))
  } else {
    counts <- colSums(data[,responses])
  }
  if(is.null(sortedTypes)) {
    sortedCounts <- sort(counts,decreasing=TRUE) ## most frequent types first  
    sortedTypes <- names(sortedCounts)
  }
  fullsortedTypes <- sortedTypes
  sure_depth <- min(c(length(sortedTypes)-1,depth))
  seq_sure_depth <- seq_len(sure_depth)
  if (input=="types") { ## contents are types
    resu <- lapply(seq_sure_depth, function(v){
      npos <- apply(data[,responses]==sortedTypes[1],1,sum)
      nneg <- apply(data[,responses],c(1,2),function(v) {v  %in% sortedTypes[-1]})
      nneg <- apply(nneg,1,sum)    
      locdata <- cbind(npos,nneg,data) 
      names(locdata)[1:2] <- binResponse
      sortedTypes <<- sortedTypes[-1]
      locdata[npos+nneg>0,] ## discards "missing data"
    })
  } else { ## colnames are types
    resu <- lapply(seq_len(min(c(length(sortedTypes)-1,depth))), function(v){
      npos <- data[,sortedTypes[1],drop=FALSE]
      nneg <- data[,sortedTypes[-1],drop=FALSE]
      nneg <- apply(nneg,1,sum)    
      locdata <- cbind(npos,nneg,data) 
      names(locdata)[1:2] <- binResponse
      sortedTypes <<- sortedTypes[-1]
      locdata[npos+nneg>0,]
    })    
  }
  names(resu) <- fullsortedTypes[seq_sure_depth] ## the names of the elements of the list; cosmetic ?
  othername <- setdiff(make.unique(c(names(resu),"others")),names(resu)) ## simply "others" unless it is already in names(resu)
  attr(resu,"sortedTypes") <- c(names(resu),othername) 
  attr(resu,"responses") <- responses
  resu ## a list of data.frames
}
