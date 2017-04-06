fit_as_ZX <- function(X.pv, # only for ncol() and colnames()... 
                      ZAL, 
                      y, ## could be taken fom processed ? 
                      n_u_h, H_global_scale, lambda_est, muetablob ,off, maxit.mean, etaFix,
                      wranefblob, processed,
                      ## supplement for ! LMM
                      eta, 
                      ## supplement for ! GLMM
                      u_h, v_h, for_init_z_args, w.resid, phi_est, 
                      ## supplement for LevenbergM
                      beta_eta,
                      ## supplement for intervals
                      for_intervals
) {
  pforpv <- ncol(X.pv)
  nobs <- length(y)
  ZAL_scaling <- 1/sqrt(wranefblob$w.ranef*H_global_scale) ## Q^{-1/2}/s
  elmts_affected_cols <- NULL
  Xscal <- make_Xscal(ZAL, ZAL_scaling=ZAL_scaling, AUGI0_ZX=processed$AUGI0_ZX,n_u_h=n_u_h)
  if (inherits(Xscal,"Matrix")) {
    ## def_sXaug_Eigen_sparse_QR calls lmwith_sparse_QRp(,pivot=FALSE)
    ## def_sXaug_Eigen_sparse_QRP calls lmwith_sparse_QRp(,pivot=TRUE)
    mMatrix_method <- .spaMM.data$options$Matrix_method
  } else mMatrix_method <- .spaMM.data$options$matrix_method
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
    old_Vscaled_beta <- c(v_h/ZAL_scaling ,beta_eta)
  } else if ( ! is.null(for_intervals)) Vscaled_beta <- c(v_h/ZAL_scaling ,for_intervals$beta_eta)
  if ( ! LMMbool) {
    constant_zAug_args <- list(n_u_h=n_u_h, nobs=nobs, pforpv=pforpv, y=y, off=off, ZAL=ZAL, processed=processed)
    if ( ! GLMMbool) {
      constant_init_z_args <- c(list(lcrandfamfam=lcrandfamfam, nobs=nobs, lambda_est=lambda_est, ZAL=ZAL),  
                                # fit_as_ZX args specific for ! GLMM:
                                for_init_z_args,
                                #
                                processed[c("cum_n_u_h","rand.families","stop.on.error")])
      constant_u_h_v_h_args <- c(processed[c("rand.families","cum_n_u_h")],
                                 processed$u_h_info, ## elements of u_h_info as elements of constant_u_h_v_h_args  
                                 list(lcrandfamfam=lcrandfamfam))
      updateW_ranefS_arglist <- c(processed[c("cum_n_u_h","rand.families")],list(lambda=lambda_est))
    } 
  } 
  ################ L O O P ##############
  for (innerj in 1:maxit.mean) {
    weight_X <- sqrt(H_global_scale*w.resid) ## sqrt(s^2 W.resid)
    ## Xscal varies within loop if !GLMM
    sXaug <- do.call(mMatrix_method,
                     list(Xaug=Xscal, weight_X=weight_X, w.ranef=wranefblob$w.ranef, H_global_scale=H_global_scale))
    if (LMMbool) {
      zAug <- c(rep(0,n_u_h), y-off)
      wzAug <- c(1/ZAL_scaling,sqrt(w.resid*H_global_scale)) *zAug ## *not* using weight_X 
    } else {
      if (LevenbergM) {
        newlik <- calc_APHLs_from_ZX(processed=processed, which=p_v_obj, sXaug=sXaug, phi_est=phi_est, 
                                         lambda_est=lambda_est, dvdu=wranefblob$dvdu, u_h=u_h, mu=muetablob$mu)[[p_v_obj]]
        if (innerj>1L) {
          gainratio <- (newlik!=-Inf) ## -Inf occurred in binary probit with extreme eta... 
          if (gainratio) { 
            summand <- LevMarblob$dVscaled_beta*(LevMarblob$rhs+ LevMarblob$dampDpD * LevMarblob$dVscaled_beta) 
            ## The two terms of the summand should be positive. In part. conv_dbetaV*rhs should be positive. 
            ## However, numerical error may lead to <0 or even -Inf
            ## Further, if there are both -Inf and +Inf elements the sum is NaN.
            summand[summand<0] <- 0
            denomGainratio <- sum(summand)
            gainratio <- 2*(newlik-oldlik)/denomGainratio ## cf computation in MadsenNT04 below 3.14, but generalized for D' != I ; numerator is reversed for maximization
          } 
          if (gainratio > 0) { 
            oldlik <- newlik
            old_Vscaled_beta <- Vscaled_beta
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
        if ( ! GLMMbool) {
          # arguments for init_resp_z_corrections_new called in calc_zAug_not_LMM
          init_z_args <- c(constant_init_z_args,
                           list(w.ranef=wranefblob$w.ranef, u_h=u_h, v_h=v_h, dvdu=wranefblob$dvdu, 
                                sXaug=sXaug, w.resid=w.resid))
        } else init_z_args <- NULL
        calc_zAug_args <- c(constant_zAug_args,
                            list(eta=eta, muetablob=muetablob, dlogWran_dv_h=wranefblob$dlogWran_dv_h, 
                                 sXaug=sXaug, ## .get_qr may be called for Pdiag calculation
                                 w.ranef=wranefblob$w.ranef, 
                                 w.resid=w.resid,
                                 init_z_args=init_z_args) )
        zAug <- do.call(".calc_zAug_not_LMM",calc_zAug_args)
        wzAug <- c(1/ZAL_scaling,sqrt(w.resid*H_global_scale)) *zAug ## *not* using weight_X 
      }
    }
    ## keep name 'w'zAug to emphasize the distinct weightings  of zaug and Xaug (should have been so everywhere)
    if ( ! is.null(for_intervals)) {
      loc_logLik_args <- list(sXaug=sXaug, processed=processed, phi_est=phi_est,
                              lambda_est=lambda_est, dvdu=wranefblob$dvdu, u_h=u_h, mu=muetablob$mu)
      newlik <- unlist(do.call("calc_APHLs_from_ZX",loc_logLik_args)[for_intervals$likfn]) # unlist keeps name
      currentDy <- (for_intervals$fitlik-newlik)
      if (currentDy <0) .warn_intervalStep(newlik,for_intervals)
      intervalBlob <- .intervalStep_ZX(old_Vscaled_beta=Vscaled_beta,
                                       sXaug=sXaug,szAug=wzAug,
                                       for_intervals=for_intervals,
                                       currentlik=newlik,currentDy=currentDy)
      Vscaled_beta <- intervalBlob$Vscaled_beta
    } else if (maxit.mean>1L && LevenbergM) { 
      oldsXaug <- sXaug
      LM_wzAug <- wzAug - as.vector(sXaug %*% old_Vscaled_beta)
      LevMarblob <- get_from_MME(sXaug=sXaug, which="LevMar_step", LMrhs=.crossprod(sXaug, LM_wzAug), damping=damping)
      Vscaled_beta <- old_Vscaled_beta + LevMarblob$dVscaled_beta 
    } else {
      Vscaled_beta <- get_from_MME(sXaug,szAug=wzAug) # vscaled= v scaling so that v has 'scale' H_global_scale
    }
    # drop, not as.vector(): names are then those of (final) eta and mu -> used by predict() when no new data
    fitted <- drop(Xscal %*% Vscaled_beta) ## length nobs+nr ! 
    fitted[ypos] <- eta <- fitted[ypos] + off
    if (is.null(etaFix$v_h)) {
      v_h <- Vscaled_beta[seq_n_u_h] * ZAL_scaling
      if (GLMMbool) {
        u_h <- v_h ## keep input wranefblob since lambda_est not changed
      } else {
        u_h_v_h_from_v_h_args <- c(constant_u_h_v_h_args,list(v_h=v_h))
        u_h <- do.call(".u_h_v_h_from_v_h",u_h_v_h_from_v_h_args)
        if ( ! is.null(attr(u_h,"v_h"))) { ## second test = if u_h_info$upper.v_h or $lower.v_h non NULL
          v_h <- attr(u_h,"v_h")
        }
        ## update functions u_h,v_h
        wranefblob <- do.call(".updateW_ranefS",c(updateW_ranefS_arglist,list(u_h=u_h,v_h=v_h)))
        ZAL_scaling <- 1/sqrt(wranefblob$w.ranef*H_global_scale) ## Q^{-1/2}/s
        if (inherits(ZAL,"Matrix")) {
          scaledZAL <- Matrix_times_Dvec(ZAL,ZAL_scaling) 
          if (is.null(elmts_affected_cols)) {
            elmts_affected_cols <- seq_len(Xscal@p[n_u_h+1L])
            which_i_affected_rows <- which(Xscal@i[elmts_affected_cols]>(n_u_h-1L))
          }
          Xscal@x[which_i_affected_rows] <- scaledZAL@x
          #Xscal[n_u_h+seq(nobs),seq_n_u_h] <- scaledZAL
        } else {
          scaledZAL <- ZAL %*% diag(x=ZAL_scaling)
          Xscal[n_u_h+seq(nobs),seq_n_u_h] <- scaledZAL
        }
      }
    } else v_h <- etaFix$v_h
    muetablob <- muetafn(eta=eta,BinomialDen=processed$BinomialDen,processed=processed) 
    if ( ! LMMbool ) w.resid <- calc.w.resid(muetablob$GLMweights,phi_est)
    beta_eta <- Vscaled_beta[n_u_h+seq_len(pforpv)]
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
  return(list(sXaug=sXaug, 
              ## used by calc_APHLs_from_ZX: (in particular can use Vscaled values contained in fitted)
              fitted=fitted, zAug=zAug, weight_X=weight_X, nobs=nobs, pforpv=pforpv, seq_n_u_h=seq_n_u_h, u_h=u_h, 
              muetablob=muetablob, 
              lambda_est=lambda_est,
              phi_est=phi_est,
              ## used by other code
              beta_eta=beta_eta, w.resid=w.resid, wranefblob=wranefblob, 
              v_h=v_h, eta=eta, innerj=innerj))
}

.intervalStep_ZX <- function(old_Vscaled_beta,sXaug,szAug,currentlik,for_intervals,currentDy) {
  #print((control.HLfit$intervalInfo$fitlik-currentlik)/(control.HLfit$intervalInfo$MLparm-old_betaV[parmcol]))
  ## voir code avant 18/10/2014 pour une implem rustique de VenzonM pour debugage  
  ## somewhat more robust algo (FR->FR: still improvable ?), updates according to a quadratic form of lik near max
  ## then target.dX = (current.dX)*sqrt(target.dY/current.dY) where dX,dY are relative to the ML x and y 
  ## A nice thing of this conception is that if the target lik cannot be approached, 
  ##   the inferred x converges to the ML x => this x won't be recognized as a CI bound (important for locoptim) 
  parmcol_ZX <- for_intervals$parmcol_ZX
  Vscaled_beta <- rep(NA,length(old_Vscaled_beta))
  if (currentDy <0) { 
    Vscaled_beta[parmcol_ZX] <- old_Vscaled_beta[parmcol_ZX]
  } else {
    currentDx <- (old_Vscaled_beta[parmcol_ZX]-for_intervals$MLparm)
    targetDy <- (for_intervals$fitlik-for_intervals$targetlik)
    Dx <- currentDx*sqrt(targetDy/currentDy)
    ## pb is if Dx=0 , Dx'=0... and Dx=0 can occur while p_v is still far from the target, because other params have not converged.
    ## FR->FR patch:
    if (currentDy<targetDy) { ## we are close to the ML: we extrapolate a bit more confidently
      min_abs_Dx <- for_intervals$asympto_abs_Dparm/1000
    } else min_abs_Dx <- 1e-6 ## far from ML: more cautious move our of Dx=0
    Dx <- sign(currentDx)*max(abs(Dx),min_abs_Dx)
    Vscaled_beta[parmcol_ZX] <- for_intervals$MLparm + Dx 
  }
  locsXaug <- sXaug[,-(parmcol_ZX),drop=FALSE]
  locszAug <- as.matrix(szAug-sXaug[,parmcol_ZX]*Vscaled_beta[parmcol_ZX])
  Vscaled_beta[-(parmcol_ZX)] <- get_from_MME(locsXaug,szAug=locszAug) 
  return(list(Vscaled_beta=Vscaled_beta)) # levQ ispresumably always dense
}
