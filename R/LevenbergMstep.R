LevenbergMstep <-
function(wAugX,LM_wAugz,damping,stop.on.error=TRUE) {
  rhs <- crossprod(wAugX,LM_wAugz) ## t(wAugX) %*% LM_wAugz ## backward: A'.z* ## rhs needed further below to assess gainratio
  ApAdDpD <- crossprodCpp(wAugX) ## t(wAugX) %*% wAugX ## A'.A=R'.R 
  dampDpD <- damping * diag(ApAdDpD) ## faster than computing qr.R(SQR ) and using it specially for this computation... ## vector
  diag(ApAdDpD) <- diag(ApAdDpD) + dampDpD  ## A'.A + damping*diag(A'.A)
  ## ## code that oddly solved B y = rhs= A'z* for B:=(A'.A + damping*diag(A'.A)) by solving B'B y =B'(A'z*) is in version 040213. For damping=0, this solved A.A'.A'.A y= A.A'.A'z*...
  ### attempt to solve 'normal equations' directly ## LMcorrected useful only here or for the ginv() call
  dbetaV <- try(solve(ApAdDpD,rhs),silent=TRUE)
  if (class(dbetaV)=="try-error") {
    ### then QR, using a standard trick: the solution of (A'A+damping D'D)y = -A'f is the solution of extended system (A // sqrt(damping) D) y = -( f // 0)
    corrD <- sqrt(dampDpD) ##  D'.D = damping* diag(A'.A) (computing diag A'A through qrR(A) is slower) ## vector
    trick <- rbind(wAugX,diag(corrD)) ## matrix
    dbetaV <- safesolve.qr.vector(qr(trick),c(LM_wAugz,rep(0,length(corrD))),stop.on.error=stop.on.error) ## 
    ## and there is nmore to this as explained by Mor'e and illustrated by the following code
    #set.seed(123)
    #m <- matrix(runif(50),ncol=5)
    #J <- Matern.corr(m,nu=1);qrJ<- qr(J);qrJ$pivot
    #b <-rnorm(10)
    #H <- t(J) %*% J
    #D <- sqrt(diag(H))
    #qrP <- qr(rbind(qr.R(qrJ),diag(sqrt(0.001)*D))) ## that's the RHS of (3.4)
    #tW <- qr.Q(qrP) # 3.5 de MOr\'e means tW.(R_lam //0) = (R//D_lam)
    #u <- t(tW) %*% c(t(qr.Q(qrJ)) %*% b,rep(0,5)) ## Mor\'e's Q ist(Q)
    #solve(qr.R(qrP),u[1:5])
    # = :
    #solve(H+0.001*diag(D^2) ,t(J) %*% b)
    ## but it misses the permutations to be generally valid !!!!!!!!!!!!!!!!!!!!!!          
  }
  if (class(dbetaV)=="try-error") {
    dbetaV <- ginv(ApAdDpD) %*% rhs
  }
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
