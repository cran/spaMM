.LevenbergMsolve_Matrix <- function(wAugX, LM_wAugz, damping){ 
  rhs <- as.vector(.crossprod(wAugX,LM_wAugz)) ## t(wAugX) %*% LM_wAugz ## backward: A'.z* ## rhs needed further below to assess
  AtAdDpD <- .crossprod(wAugX)
  dampDpD <- damping * diag(AtAdDpD)
  nc <- ncol(wAugX)
  diagPos <- seq.int(1L,nc^2,nc+1L)
  AtAdDpD[diagPos] <- AtAdDpD[diagPos] + dampDpD
  dbetaV <- as.vector(solve(AtAdDpD,rhs))
  return(list(dbetaV=dbetaV, dampDpD=dampDpD, rhs=rhs)) #,damping=damping)) 
}


.LevenbergMstepCallingCpp <- function(wAugX,LM_wAugz,damping) {
  ## FIXME perhaps http://eigen.tuxfamily.org/dox/unsupported/LMonestep_8h_source.html could be useful ???
  rhs <- as.vector(.crossprodCpp_d(wAugX,LM_wAugz)) ## t(wAugX) %*% LM_wAugz ## backward: A'.z* ## rhs needed further below to assess gainratio
  resu <- .LevenbergMsolveCpp(wAugX, rhs, damping ) 
  #resu <- .e_LevenbergMsolveCpp(wAugX, LM_wAugz, damping ) ## essai 07/2016, using Eigen QR, not convincing
  resu$rhs <- rhs
  #resu$damping <- damping
  return(resu) ## format: list(dbetaV, dampDpD, rhs)
}  

