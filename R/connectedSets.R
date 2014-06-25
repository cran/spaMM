connectedSets <-
function(indices) {
  uI <- unique(as.vector(indices)) ## vertex indices
  rownames(indices) <- seq_len(nrow(indices))
  simplicesTable <- list()
  it <- 1
  while (nrow(indices)>1) {
    oldvertexSet <- numeric(0)
    growingvertexSet <- uI[1] # a vertex
    while(length(oldvertexSet)!=length(growingvertexSet)) {
      oldvertexSet <- growingvertexSet
      rowIndices <- which(apply(indices,1,function(v) { any(v %in% growingvertexSet)}))
      connectedRows <- indices[rowIndices,]
      growingvertexSet <- unique(as.vector(connectedRows))
    }
    simplicesTable[[it]] <- growingvertexSet
    uI <- setdiff(uI,growingvertexSet)
    indices <- indices[-rowIndices,,drop=FALSE]
    it <- it+1
  }
  if (nrow(indices)==1L) simplicesTable[[it]] <- as.numeric(indices) ## isolated simplex
  lapply(simplicesTable,as.numeric) 
}
