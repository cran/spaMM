rvolTriangulation <-
function(n=1,volTriangulationObj,replace=TRUE) { 
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
}
