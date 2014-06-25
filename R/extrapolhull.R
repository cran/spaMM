extrapolhull <-
function(Vhull,extrapol=1.4) { ##
  if (nrow(Vhull) <2) return(matrix(nrow=0,ncol=ncol(Vhull)))
  ## ELSE
  vertexbary <- colMeans(Vhull) ## not bary for uniform density, except for simplices
  ## extremely lightweight solution: random extrapol along all directions defined by the vertices and the bary
  extrap <- runif(nrow(Vhull))*(extrapol-1)+1
  t((t(Vhull)-vertexbary)* extrap+ vertexbary)  ## bary+ extrapol*(v-bary)
}
