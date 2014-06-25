volTriangulation <-
function(vertices) { ## by contrast to barycenter of vertices
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
