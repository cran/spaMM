sampleGridFromLowUp <-
function(LowUp,n=NULL,gridSteps=NULL,sampling=NULL) {
  ## grid from lower, upper
  lower <- LowUp$lower
  upper <- LowUp$upper
  d <- length(lower)
  if (is.null(sampling)) {
    if (d<4) {sampling="grid"} else {sampling="rvolTriangulation"}
  }
  byvar <- t(rbind(unlist(lower),unlist(upper))) 
  byvar <- 0.999 * byvar + 0.001 *rowMeans(byvar)
  grillelist <- list()
  if (is.null(gridSteps)) {
    gridSteps <- c(20,6,4,3)[min(d,4)] ## => 20  36  64  81 243 729 points 
  }
  for(name in rownames(byvar)) {grillelist[[name]] <- seq(byvar[name,1],byvar[name,2],length.out=gridSteps)}
  pargrid <- expand.grid(grillelist)
  if (sampling=="rvolTriangulation") { ## while randomly sample the simplices defined by the regular grid
    vT <- volTriangulation(pargrid) ## note no convhulln call -> many internal simplices
    if (is.null(n)) n <- gridSteps^min(d,6) ## default maxi 729 for high d
    pargrid <- rbind(pargrid,rvolTriangulation(n,vT)) ## regular + random
  } else { ## more structured sampling: will sample the grid 
    attr(pargrid,"regularGrid") <- seq_len(nrow(pargrid))
    ## redefines a grid
    insides <- lapply(grillelist, function(v) {(v[2]-v[1])/2+v[-c(gridSteps)]}) ## 19 5 3 2 2 2...
    stepsizes <- unlist(lapply(insides, function(v) {v[2]-v[1]})) 
    insides <- expand.grid(insides) ## 19 25 27 16 32 64 ...
    randgrid <- sampleNearby(insides,n=nrow(insides),stepsizes=stepsizes) ## not really nearby given the large stepsizes
    pargrid <- rbind(pargrid,randgrid) ## regular + random
  }
  ##
  pargrid
}
