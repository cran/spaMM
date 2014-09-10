sampleNextPoints <-
function(n=1,optr,minPtNbr=1,##sqrt(nrow(Xpredy)+1),
                             D.resp=NULL,replace=TRUE,expand=1) { ## different conception, more adaptive
  Xpredy <- optr$predictions
  fittedPars <- attr(Xpredy,"fittedPars")
  vT <- volTriangulation(as.matrix(Xpredy[,fittedPars])) 
  if (is.infinite(expand)) {
    upperSimplices <- seq_len(nrow(vT$simplicesTable))
    innerVertexIndices <- seq_len(nrow(Xpredy)) 
  } else {
    ## seek simplices which involve 'interesting' points
    innerVertexIndices <- upperPoints(Xpredy,minPtNbr=minPtNbr,D.resp=D.resp) ## indices; FR->FR maybe not most efficient as Xpredy is already sorted
    ## the following defines a NON CONVEX set
    upperSimplices <- apply(vT$simplicesTable,1,function(v) { any(v %in% innerVertexIndices)}) ## indices !
    for (it in seq_len(expand-1)) {
      innerVertexIndices <- unique(as.vector(vT$simplicesTable[upperSimplices,]))
      upperSimplices <- apply(vT$simplicesTable,1,function(v) { any(v %in% innerVertexIndices)}) ## indices !    
    }
  }
  subvT <- vT
  subvT$innerVertexIndices <- innerVertexIndices
  subvT$simplicesTable <- vT$simplicesTable[upperSimplices,]
  subvT$vol <- vT$vol[upperSimplices]
  ## but keep the original $vertices since the subvT$simplicesTable still refer to original rows 
  resu <- rvolTriangulation(n=n,subvT,replace=replace)
  colnames(resu) <- fittedPars
  resu <- data.frame(resu)
  attr(resu,"info") <- subvT
  return(resu)
}
