LevenbergMstepCallingCpp <- function(wAugX,LM_wAugz,damping) {
  ## FR->FR perhaps http://eigen.tuxfamily.org/dox/unsupported/LMonestep_8h_source.html could be useful ???
  rhs <- drop(crossprod(wAugX,LM_wAugz)) ## t(wAugX) %*% LM_wAugz ## backward: A'.z* ## rhs needed further below to assess gainratio
  if (TRUE) {
    resu <- LevenbergMsolveCpp(wAugX, rhs, damping ) 
  } else resu <- e_LevenbergMsolveCpp(wAugX, LM_wAugz, damping ) ## essai 07/2016
  resu$rhs <- rhs
  return(resu) 
}  

# R version of e_LevenbergMsolveCpp: algo using QR of extended system (cf Moré or Nocedal & Wright):
e_LevenbergMstep <- function(wAugX,LM_wAugz,damping) {
  rhs <- drop(crossprod(wAugX,LM_wAugz)) ## t(wAugX) %*% LM_wAugz ## backward: A'.z* ## rhs needed further below to assess gainratio
  QR_A <- Rcpp_QR(wAugX) ## should be precomputed once to save time and given as input to the function
  dampDpD <- damping * colSums(wAugX*wAugX)
  corrD <- sqrt(dampDpD)
  QR_QR_R <- Rcpp_QR(rbind(QR_A$R,diag(corrD,nrow=length(corrD)))) ## a optimiser selon HOuseholder et Givens ?
  nc <- ncol(wAugX)
  z <- rbind(t(QR_A$Q) %*% LM_wAugz,matrix(rep(0,nc),ncol=1))
  z <- t(QR_QR_R$Q[,1:nc]) %*% z ## better to get directy the thin Q ?
  dbetaV <- backsolve(QR_QR_R$R[1:nc,],z)
  if (inherits(dbetaV,"try-error")) {
    dbetaV <- ginv(QR_QR_R$R[1:nc,]) %*% z
  }
  resu <- list(dbetaV=dbetaV,dampDpD=dampDpD,rhs=rhs)
  return(resu) 
}  



## here version 1.5.3 has an interesting signed.wAugX concept
LevenbergMstep <- function(wAugX,LM_wAugz,damping,stop.on.error=TRUE) {
  rhs <- drop(crossprod(wAugX,LM_wAugz)) ## t(wAugX) %*% LM_wAugz ## backward: A'.z* ## rhs needed further below to assess gainratio
  if (inherits(wAugX,"Matrix")) {
    ApAdDpD <- crossprod(wAugX)
  } else ApAdDpD <- crossprodCpp(wAugX) ## t(wAugX) %*% wAugX ## A'.A=R'.R     
  dampDpD <- damping * diag(ApAdDpD) ## faster than computing qr.R(qr_wAugX ) and using it specially for this computation... ## vector
  nc <- ncol(ApAdDpD)
  diagPos <- seq.int(1L,nc^2,nc+1L)
  ApAdDpD[diagPos] <- ApAdDpD[diagPos] + dampDpD ## A'.A + damping*diag(A'.A)
  ### attempt to solve 'normal equations' directly ## LMcorrected useful only here or for the ginv() call
  dbetaV <- try(solve(ApAdDpD,rhs),silent=TRUE)
  if (inherits(dbetaV,"try-error")) {
    ### then QR, using a standard trick: the solution of (A'A+damping D'D)y = -A'f is the solution of extended system (A // sqrt(damping) D) y = -( f // 0)
    corrD <- sqrt(dampDpD) ##  D'.D = damping* diag(A'.A) (computing diag A'A through qrR(A) is slower) ## vector
    trick <- rbind(wAugX,diag(corrD)) ## matrix
    dbetaV <- safesolve.qr.vector(qr(trick),c(LM_wAugz,rep(0,length(corrD))),stop.on.error=stop.on.error) ## 
    ## see comments on More77 suggesting that the QR of the perturbed matrix involves the R of the unperturbed one 
  }
  if (inherits(dbetaV,"try-error")) {
    dbetaV <- ginv(ApAdDpD) %*% rhs
  }
  dbetaV <- as.numeric(dbetaV) ## may be dgeMatrix if wAugX was dgCMatrix
  return(list(dbetaV=dbetaV,rhs=rhs,dampDpD=dampDpD))
}

