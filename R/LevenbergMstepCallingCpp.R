LevenbergMstepCallingCpp <-
function(wAugX,LM_wAugz,damping) {
  ## FR->FR perhaps http://eigen.tuxfamily.org/dox/unsupported/LMonestep_8h_source.html could be useful ???
  rhs <- crossprod(wAugX,LM_wAugz) ## t(wAugX) %*% LM_wAugz ## backward: A'.z* ## rhs needed further below to assess gainratio
  resu <- LevenbergMsolveCpp(wAugX, rhs, damping )
  dbetaV <- as.matrix(resu$dbetaV)
  dampDpD <- resu$dampDpD
  summand <- dbetaV*(rhs+ dampDpD * dbetaV) 
  ## In the summand, all terms should be positive. dbetaV*rhs should be positive. However, numerical error may lead to <0 or even -Inf
  ## Further, if there are both -Inf and +Inf elements the sum is NaN and HLfit fails.
  summand[summand<0] <- 0
  denomGainratio <- sum(summand)
  #if (is.nan(denomGainratio)) {
  #  mess <- pastefrom("NaN 'denomGainratio'.",prefix="(!) From ") 
  #  stop(mess)
  #}
  return(list(dbetaV=dbetaV,denomGainratio=denomGainratio)) 
}
