rsimplex <-
function(simplex) { ## to sample uniformly wrt volume within a simplex, not uniformly wrt perimeter
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
