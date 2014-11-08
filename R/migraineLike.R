#load("../predictions.Rdata")
#attr(predictions,"fittedPars") <- colnames(predictions)[1:2]
#attr(predictions,"respName") <- colnames(predictions)[3]
#attr(predictions,"RMSy") <- 2

#require(geometry) ## delaunayn,convhulln


upperPoints <- function(Xpredy,minPtNbr=1#sqrt(nrow(Xpredy)+1)
                        ,D.resp) {
  respName <- attr(Xpredy,"respName")
  orderResp <- order(Xpredy[,respName], decreasing = TRUE)
  orderedXy <- Xpredy[orderResp,]
  Dy <- orderedXy[1,respName]-orderedXy[,respName]
  ptNbr <- max(c(which(Dy<D.resp)),minPtNbr) 
  orderResp[seq_len(ptNbr)] ## now returns indices
}

#nH <- upperPoints(predictions)

elim.redundant.V <- function(vertices) { ## removes redundant vertices
  if (nrow(vertices)<=ncol(vertices)) { ## convhulln crashes!
    minimalvertices<-vertices ## stupid case, should never occur
  } else if (ncol(vertices)==1) {
    minimalvertices<-array(c(min(vertices),max(vertices)),dim=c(2,1))
  } else {
    minimalvertices<-vertices[unique(as.numeric(convhulln(vertices, "Pp"))),] ## removes redundant vertices
  }
  colnames(minimalvertices)<-colnames(vertices)
  minimalvertices
}

#mvH <- elim.redundant.V(nH)

volTriangulation <- function(vertices) { ## by contrast to barycenter of vertices
  tc <- delaunayn(vertices,"Pp") ## triangulation by simplices
  pmul <- cbind(-1,diag(rep(1,ncol(vertices))))
  vb <- apply(tc,1,function(v){
    simplex <- vertices[v,,drop=FALSE]
    # pmul %*% simplex is simplex[-1,]-simplex[1,] for each simplex
    # abs(det()) is its volume up to a scaling factor constant among simplices
    vol <- abs(det(pmul %*% simplex))
    bary <- colMeans(simplex)
    c(vol,bary) ## volume, and barycenter \n\
  }) 
  resu <- list(vol=vb[1,],vertices=vertices,simplicesTable=tc)
  class(resu) <- c("volTriangulation",class(resu))
  resu
}

#vT <- volTriangulation(mvH)

rsimplex <-function(simplex) { ## to sample uniformly wrt volume within a simplex, not uniformly wrt perimeter
  ## even extrapolated values have a high chance of falling into an adjacent simplex
  bary <- colMeans(simplex) ## OK pour bary de d-dim simplex avec densité uniforme 
  #simplex <- t(extrapol*(t(simplex)-bary)+ bary)  ## bary+ extrapol*(v-bary)
  d <- NROW(simplex)-1 ## 2 pour triangle
  if(NCOL(simplex)!= d) {stop("(!) From 'rsimplex': simplex should have one more row than columns.")}
  u <- runif(d)
  upow <- u^(1/seq(d,1)) ## u_1^{1/2},u_2
  ## loop eliminates one dimension at each step:
  for(it in seq_len(d)) simplex <- t(upow[it]*(t(simplex[-1,,drop=FALSE])-simplex[1,])+ simplex[1,])  ## v_1+ upow[it](v[-1]-v_1)
  simplex ## a point
}


## ideally we would have a rvolume S3 generic with .volTriangulation and .data.frame methods
rvolTriangulation <- function(n=1,volTriangulationObj,replace=TRUE) { 
  if (! inherits(volTriangulationObj,"volTriangulation")) stop("(!) From 'volTriangulation': input 'volTriangulationObj' is not a volTriangulation object.")
  simplexProbs <- volTriangulationObj$vol
  simplexProbs <- simplexProbs/sum(simplexProbs)
  whichSimplex <- sample(seq_len(length(simplexProbs)),n,prob=simplexProbs,replace=replace)
  vertices <- volTriangulationObj$vertices
  vI <- volTriangulationObj$simplicesTable
  resu <- sapply(whichSimplex,function(idx) {rsimplex(vertices[vI[idx,],,drop=FALSE])})
  if (ncol(vertices)==1L) {
    resu <- as.matrix(resu)
  } else resu <- t(resu)
#  colnames(resu) <- colnames(vertices) ## currently no name info in the arguments of rvolTriangulation!
  return(resu)
} ## end def rvolDelaunayn

#rvolTriangulation(10,vT)
#plot(vT$vertices)
#points(rvolTriangulation(1000,vT),col="red")


old.nextPoints <- function(n=1,optr,replace=TRUE) { ## random sampling of volume defined from previous fit
  uP <- upperPoints(optr$predictions) ## indices
  uP <- optr$predictions[uP,attr(optr$predictions,"fittedPars")]
  uP <- rbind(uP,optr$par) ## not sure this is useful for volumetric sampling
  erV <- elim.redundant.V(uP)  
  vT <- volTriangulation(erV)
  rvolTriangulation(n,vT,replace=replace)
}


# Retruns sampled locations
# attribute "info" contains 
#    the original $vertices (indices), 
#    $upperVertexIndices the indices of the best points, 
#    $simplicesTable the subset of simplices involving the best points
#    $vol the volumes of subset of simplices involving the best points
# default minPtNbr affect exploration of (provisionally) suboptimal peaks
# expand=2 allows expansion at least towards a local peak
sampleNextPoints <- function(n=1,optr,minPtNbr=1,##sqrt(nrow(Xpredy)+1),
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

connectedSets <- function(indices) {
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
} ## returns point indices

extrapolhull <-function(Vhull,extrapol=1.4) { ##
  if (nrow(Vhull) <2) return(matrix(nrow=0,ncol=ncol(Vhull)))
  ## ELSE
  vertexbary <- colMeans(Vhull) ## not bary for uniform density, except for simplices
  ## extremely lightweight solution: random extrapol along all directions defined by the vertices and the bary
  extrap <- runif(nrow(Vhull))*(extrapol-1)+1
  t((t(Vhull)-vertexbary)* extrap+ vertexbary)  ## bary+ extrapol*(v-bary)
}  