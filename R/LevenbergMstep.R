.LevenbergMsolve_Matrix <- function(wAugX, LM_wAugz, damping){ 
  rhs <- as.vector(Matrix::crossprod(wAugX,LM_wAugz)) ## t(wAugX) %*% LM_wAugz ## backward: A'.z* ## rhs needed further below to assess
  AtAdDpD <- Matrix::crossprod(wAugX)
  dampDpD <- damping * diag(AtAdDpD)
  nc <- ncol(wAugX)
  diagPos <- seq.int(1L,nc^2,nc+1L)
  AtAdDpD[diagPos] <- AtAdDpD[diagPos] + dampDpD
  dbetaV <- as.vector(solve(AtAdDpD,rhs))
  return(list(dbetaV=dbetaV, dampDpD=dampDpD, rhs=rhs)) #,damping=damping)) 
}


.LevenbergMstepCallingCpp <- function(wAugX,LM_wAugz,damping) {
  ## FR->FR perhaps http://eigen.tuxfamily.org/dox/unsupported/LMonestep_8h_source.html could be useful ???
  rhs <- as.vector(.crossprodCpp(wAugX,LM_wAugz)) ## t(wAugX) %*% LM_wAugz ## backward: A'.z* ## rhs needed further below to assess gainratio
  if (TRUE) {
    resu <- .LevenbergMsolveCpp(wAugX, rhs, damping ) 
  } else resu <- .e_LevenbergMsolveCpp(wAugX, LM_wAugz, damping ) ## essai 07/2016
  resu$rhs <- rhs
  resu$damping <- damping
  return(resu) 
}  

