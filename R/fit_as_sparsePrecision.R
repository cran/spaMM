fit_as_sparsePrecision <- function(X.pv, 
                                   ZAL, ## avoid ## still used by calc_zAug_not_LMM ## hmmm and ZALX
                                   y, n_u_h, 
                       H_global_scale=1, 
                       lambda_est, muetablob,
                       off, maxit.mean, etaFix,
                       ## for ! LMM
                       eta, 
                       ## supplement for LevenbergM
                       beta_eta,
                       ## supplement for ! GLMM
                       wranefblob, u_h, v_h, w.resid, phi_est,
                       for_init_z_args,
                       for_intervals,
                       ##
                       processed,rho) {
  
  pforpv <- ncol(X.pv)
  nobs <- length(y)
  elmts_affected_cols <- NULL
  #ZX <- cbind2(processed$AUGI0_ZX$ZAfix,X.pv)
  ZALX <- cbind2(ZAL,X.pv)
  #Xscal <- make_Xscal(ZAL, ZAL_scaling=ZAL_scaling, AUGI0_ZX=processed$AUGI0_ZX,n_u_h=n_u_h)
  seq_n_u_h <- seq_len(n_u_h)
  ypos <- n_u_h+seq_len(nobs)
  lcrandfamfam <- attr(processed$rand.families,"lcrandfamfam")
  LMMbool <- processed$LMMbool
  GLMMbool <- processed$GLMMbool
  if (maxit.mean>1L) conv.threshold <- processed$conv.threshold
  LevenbergM <- (processed$LevenbergM && is.null(for_intervals))
  if (LevenbergM) { 
    damping <- 1e-7 ## cf Madsen-Nielsen-Tingleff; high valuesgive poor results
    dampingfactor <- 2
    if (processed$HL[1L]==0L) {p_v_obj <-"hlik"} else p_v_obj <-"p_v" ## objective for beta(_v) estim only: != outer obj 
    #old_Vscaled_beta <- c(v_h/ZAL_scaling ,beta_eta)
    old_v_h_beta <- list(v_h=v_h,beta_eta=beta_eta)
    old_eta <- X.pv %*% old_v_h_beta$beta_eta + ZAL %*% old_v_h_beta$v_h
  } else if ( ! is.null(for_intervals)) v_h_beta <- list(v_h=v_h, beta_eta=for_intervals$beta_eta) #stop("Interval calculation not implemented for sparse-precision matrix methods.") # FIXME
  if ( ! LMMbool) {
    constant_zAug_args <- list(n_u_h=n_u_h, nobs=nobs, pforpv=pforpv, y=y, off=off, ZAL=ZAL, processed=processed)
    if ( ! GLMMbool) {
      stop("sparse-precision matrix methods work only for GLMMs")
    } 
  } 
  ################ L O O P ##############
  for (innerj in 1:maxit.mean) {
    sXaug <- do.call("def_AUGI0_ZX_sparsePrecision",
                     list(AUGI0_ZX=processed$AUGI0_ZX, rho=rho, lambda=lambda_est[1], 
                          w.resid=w.resid))  ## FR->FR rename ?
    if (LMMbool) {
      zInfo <- list(z2=rep(0,n_u_h),z1=y-off,sscaled=0)
    } else {
      if (LevenbergM) {
        newlik <- calc_APHLs_from_ZX(processed=processed, which=p_v_obj, sXaug=sXaug, phi_est=phi_est, 
                                     lambda_est=lambda_est, dvdu=wranefblob$dvdu, u_h=u_h, mu=muetablob$mu)[[p_v_obj]]
        if (innerj>1L) {
          gainratio <- (newlik!=-Inf) ## -Inf occurred in binary probit with extreme eta... 
          if (gainratio) { 
            dv_h_beta <- unlist(LevMarblob[c("dv_h","dbeta_eta")])
            summand <- dv_h_beta*(LevMarblob$LMrhs+ LevMarblob$dampDpD * dv_h_beta) 
            ## The two terms of the summand should be positive. In part. dv_h_beta*LMrhs should be positive. 
            ## However, numerical error may lead to <0 or even -Inf
            ## Further, if there are both -Inf and +Inf elements the sum is NaN.
            summand[summand<0] <- 0
            denomGainratio <- sum(summand)
            gainratio <- 2*(newlik-oldlik)/denomGainratio ## cf computation in MadsenNT04 below 3.14, but generalized for D' != I ; numerator is reversed for maximization
          } 
          if (gainratio > 0) { 
            oldlik <- newlik
            old_eta <- eta
            old_v_h_beta <- v_h_beta
            old_relV_beta <- relV_beta
            ## cf Madsen-Nielsen-Tingleff again, and as in levmar library by Lourakis
            damping <- damping * max(1/3,1-(2*gainratio-1)^3)  
            dampingfactor <- 2
          } else { ## unsuccesful step: do not update anything
            sXaug <- oldsXaug ## we refit with the old sXaug and an updated damping
            damping <- damping*dampingfactor
            dampingfactor <- dampingfactor*2 
          }
        } else {
          oldlik <- newlik
        }
      }
      if ( innerj==1L || ! LevenbergM || gainratio>0L ) {
        init_z_args <- NULL ## NULL only for GLMMs...
        calc_zAug_args <- c(constant_zAug_args,
                            list(eta=eta, muetablob=muetablob, dlogWran_dv_h=wranefblob$dlogWran_dv_h, 
                                 sXaug=sXaug, ## .get_qr may be called for Pdiag calculation
                                 w.ranef=wranefblob$w.ranef, 
                                 w.resid=w.resid,
                                 init_z_args=init_z_args) )
        augz <- do.call(".calc_zAug_not_LMM",calc_zAug_args)
        zInfo <- attributes(augz) # .AUGI0_ZX_sparsePrecision use only the attributes; 
        # fit_as_ZX -> other get_from_MME methods uses augz: FR->FR: modify fit_as_ZX ?
      }
    }
    ## keep name 'w'zAug to emphasize the distinct weightings  of zaug and Xaug (should have been so everywhere)
    if ( ! is.null(for_intervals)) {
      loc_logLik_args <- list(sXaug=sXaug, processed=processed, phi_est=phi_est,
                              lambda_est=lambda_est, dvdu=wranefblob$dvdu, u_h=u_h, mu=muetablob$mu)
      newlik <- unlist(do.call("calc_APHLs_from_ZX",loc_logLik_args)[for_intervals$likfn]) # unlist keeps name
      currentDy <- (for_intervals$fitlik-newlik)
      if (currentDy <0) .warn_intervalStep(newlik,for_intervals)
      intervalBlob <- .intervalStep_I0ZX(old_v_h_beta=v_h_beta,
                                        sXaug=sXaug,zInfo=zInfo,
                                        for_intervals=for_intervals,
                                        currentlik=newlik,currentDy=currentDy)
      v_h_beta <- intervalBlob$v_h_beta
    } else if (maxit.mean>1L && LevenbergM) {
      oldsXaug <- sXaug
      zInfo$z1_eta <- zInfo$z1-old_eta
      zInfo$phi_v <- old_v_h_beta$v_h/lambda_est
      LevMarblob <- get_from_MME(sXaug=sXaug, which="LevMar_step", LM_z=zInfo, damping=damping)
      ## compute $LMrhs now bc some variables ued here may be modified before $LMrhs is used 
      ## $LMrhs is the rhs in the direct solution of the full system (by chol2inv in the dense QR case)
      ## The 'betaFirst' algo in sparsePrecision does not use it, onlythe gainratio code does  
      LevMarblob$LMrhs <- c( 
        .crossprod(ZAL, attr(sXaug,"w.resid") * (zInfo$z1-old_eta)) - old_v_h_beta$v_h/lambda_est, # Z'W(z_1-eta)-Phi v_h
        .crossprod(X.pv, attr(sXaug,"w.resid") * (augz[-seq_n_u_h]-old_eta)) # X'W(z_1-sscaled-eta)
      ) ## the RHS of 1st eq in "Details of stepwise solution" less the terms in red therein (including the red terms within z1)
      v_h_beta <- list(
        v_h=old_v_h_beta$v_h + LevMarblob$dv_h,
        beta_eta=old_v_h_beta$beta_eta + LevMarblob$dbeta_eta
      )
    } else {
      v_h_beta <- get_from_MME(sXaug,szAug=zInfo) 
    }
    # drop, not as.vector(): names are then those of (final) eta and mu -> used by predict() when no new data
    fitted <- drop(ZALX %*% unlist(v_h_beta)) ## length nobs+nr ! 
    eta <- fitted + off
    if ( is.null(etaFix$v_h)) {
      v_h <- v_h_beta$v_h
      ## u_h is not used internally but returned, modified only if v_h has been modified.
      ## Hence we update u_h ONLY when v_h is updated (= not in etaFix$v_h case, whee the input u_h is kept)
      u_h <- v_h ## GLMMbool only: u=v 
    } else v_h <- etaFix$v_h
    muetablob <- muetafn(eta=eta,BinomialDen=processed$BinomialDen,processed=processed) 
    if ( ! LMMbool ) w.resid <- calc.w.resid(muetablob$GLMweights,phi_est)
    beta_eta <- v_h_beta$beta_eta
    if (innerj<maxit.mean) {
      relV_beta <- c(v_h*sqrt(wranefblob$w.ranef),beta_eta)  ## convergence on v_h relative to sqrt(lambda), more exactly for Gaussian
      break_cond <- (innerj>1L && mean(abs(relV_beta - old_relV_beta)) < conv.threshold/10)
      if (is.na(break_cond)) {
        if (anyNA(relV_beta)) {
          stop("Numerical problem: try control.HLfit=list(LevenbergM=TRUE)")
        } else stop("Error in evaluating break condition")
      } else if (break_cond) break
      if ( ! LevenbergM || innerj==1L) old_relV_beta <- relV_beta ## otherwise this is updated conditionnally above.
    } else break
  } ################ E N D LOOP ##############
  names(beta_eta) <- colnames(X.pv)
  #weight_X <- sqrt(H_global_scale*w.resid) ## sqrt(s^2 W.resid)
  return(list(sXaug=sXaug, 
              ## used by calc_APHLs_from_ZX:
              #fitted=fitted, ## FIXME: removed so that no shortcut a la Bates in calc_APHLs_from_ZX; reimplement the shorcut in that fn?  
              weight_X=NA, nobs=nobs, pforpv=pforpv, 
              seq_n_u_h=seq_n_u_h, 
              v_h=v_h, 
              u_h=u_h,
              muetablob=muetablob, 
              lambda_est=lambda_est,
              phi_est=phi_est,
              ## used by other code
              beta_eta=beta_eta, w.resid=w.resid, wranefblob=wranefblob, 
              eta=eta, innerj=innerj))
} 

.intervalStep_I0ZX <- function(old_v_h_beta,sXaug,zInfo,currentlik,for_intervals,currentDy) {
  #print((control.HLfit$intervalInfo$fitlik-currentlik)/(control.HLfit$intervalInfo$MLparm-old_betaV[parmcol]))
  ## voir code avant 18/10/2014 pour une implem rustique de VenzonM pour debugage  
  ## somewhat more robust algo (FR->FR: still improvable ?), updates according to a quadratic form of lik near max
  ## then target.dX = (current.dX)*sqrt(target.dY/current.dY) where dX,dY are relative to the ML x and y 
  ## A nice thing of this conception is that if the target lik cannot be approached, 
  ##   the inferred x converges to the ML x => this x won't be recognized as a CI bound (important for locoptim) 
  parmcol_ZX <- for_intervals$parmcol_ZX
  v_h_beta_vec <- old_v_h_beta_vec <- unlist(old_v_h_beta)
  v_h_beta_vec[] <- NA
  if (currentDy <0) { 
    v_h_beta_vec[parmcol_ZX] <- old_v_h_beta_vec[parmcol_ZX]
  } else {
    currentDx <- (old_v_h_beta_vec[parmcol_ZX]-for_intervals$MLparm)
    targetDy <- (for_intervals$fitlik-for_intervals$targetlik)
    Dx <- currentDx*sqrt(targetDy/currentDy)
    ## pb is if Dx=0 , Dx'=0... and Dx=0 can occur while p_v is still far from the target, because other params have not converged.
    ## FR->FR patch:
    if (currentDy<targetDy) { ## we are close to the ML: we extrapolate a bit more confidently
      min_abs_Dx <- for_intervals$asympto_abs_Dparm/1000
    } else min_abs_Dx <- 1e-6 ## far from ML: more cautious move our of Dx=0
    Dx <- sign(currentDx)*max(abs(Dx),min_abs_Dx)
    v_h_beta_vec[parmcol_ZX] <- for_intervals$MLparm + Dx 
  }
  locAUGI0_ZX <- sXaug
  parmcol_X <- for_intervals$parmcol_X
  locAUGI0_ZX$X.pv <- locAUGI0_ZX$X.pv[,-(parmcol_X),drop=FALSE]
  locAUGI0_ZX$ZeroBlock <- locAUGI0_ZX$ZeroBlock[,-(parmcol_X),drop=FALSE]
  zInfo$z1 <- zInfo$z1 - sXaug$X.pv[,parmcol_X]*v_h_beta_vec[parmcol_ZX]
  v_h_beta_vec[-(parmcol_ZX)] <- unlist(get_from_MME(locAUGI0_ZX,szAug=zInfo)) 
  return(list(v_h_beta=relist(v_h_beta_vec,old_v_h_beta))) 
}