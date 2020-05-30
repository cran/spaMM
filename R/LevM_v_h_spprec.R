.solve_v_h_IRLS_spprec <- # only for LevM && is_HL_1_1
  function(X.pv, 
           ZAL, y, ## could be taken fom processed ? 
           n_u_h, 
           #H_global_scale, 
           lambda_est, muetablob=NULL, off, maxit.mean, etaFix,
           wranefblob, processed,
           ## supplement for ! LMM
           phi_est, 
           ## supplement for for LevenbergM 
           w.resid=NULL, 
           ## supplement for LevenbergM
           beta_eta,
           ## supplement for ! GLMM (??)
           u_h, v_h, for_init_z_args, 
           #
           trace=FALSE,
           stylefn=.spaMM.data$options$stylefns$vloop,
           looseness,
           LevMarblob=NULL,
           corrPars,
           dampings_env
  ) {
    pforpv <- ncol(X.pv)
    nobs <- length(y)
    seq_n_u_h <- seq_len(n_u_h)
    ypos <- n_u_h+seq_len(nobs)
    lcrandfamfam <- attr(processed$rand.families,"lcrandfamfam")
    LMMbool <- processed$LMMbool
    GLMMbool <- processed$GLMMbool
    not_moving <- FALSE
    old_relV <- NULL
    damped_WLS_blob <- NULL
    # looseness controlled by .wrap_do_damped_WLS_outer(..., looseness, ...) which may be overwritten when "strict_v|b"
    pot_tol <- processed$spaMM_tol$v_pot_tol * looseness
    d_relV_tol <- processed$spaMM_tol$d_relV_tol * looseness
    constant_zAug_args <- list(n_u_h=n_u_h, nobs=nobs, pforpv=pforpv, y=y, off=off, ZAL=ZAL, processed=processed)
    if ( ! GLMMbool) {
      constant_init_z_args <- c(list(lcrandfamfam=lcrandfamfam, nobs=nobs, lambda_est=lambda_est, ZAL=ZAL),  
                                # fit_as_ZX args specific for ! GLMM:
                                for_init_z_args,
                                #
                                mget(c("cum_n_u_h","rand.families","stop.on.error"),envir=processed))
      delayedAssign("constant_u_h_v_h_args", 
                    c(mget(c("cum_n_u_h","rand.families"),envir=processed),
                      processed$u_h_info, ## elements of u_h_info as elements of constant_u_h_v_h_args, NOT u_h_info as element of...  
                      list(lcrandfamfam=lcrandfamfam)))
      updateW_ranefS_constant_arglist <- c(mget(c("cum_n_u_h","rand.families"),envir=processed),list(lambda=lambda_est))
    } 
    
    ##### initial sXaug
    ZAL_scaling <- 1  ## TAG: scaling for spprec
    if (is.null(muetablob)) { ## NULL input eta allows NULL input muetablob
      eta  <- off + drop(processed$AUGI0_ZX$X.pv %*% beta_eta) + drop(ZAL %id*% v_h) 
      muetablob <- .muetafn(eta=eta,BinomialDen=processed$BinomialDen,processed=processed) 
    }
    ## weight_X and Xscal varies within loop if ! LMM since at least the GLMweights in w.resid change
    if ( is.null(w.resid) ) w.resid <- .calc_w_resid(muetablob$GLMweights,phi_est)
    ## needs adjMatrix and corrPars to define Qmat
    update_sXaug_constant_arglist <- list(AUGI0_ZX=processed$AUGI0_ZX, corrPars=corrPars, 
                                          cum_n_u_h=processed$cum_n_u_h #,H_global_scale=H_global_scale
                                          ) 
    #weight_X <- .calc_weight_X(w.resid, H_global_scale) ## sqrt(s^2 W.resid)
    if (trace) cat(stylefn(".")) ## hmff blue (vloop) F I X M E
    sXaug <- do.call(processed$AUGI0_ZX$envir$method, # ie, def_AUGI0_ZX_sparsePrecision
                     c(update_sXaug_constant_arglist,
                       list(w.ranef=wranefblob$w.ranef, 
                            #weight_X=weight_X, 
                            w.resid=w.resid)))
    Vscaled_beta <- list(v_h=v_h/ZAL_scaling,beta_eta=beta_eta)
    
    break_info <- list()
    ################ L O O P ##############
    for (innerj in 1:maxit.mean) {
      ##### get the lik of the current state
      if (is.null(damped_WLS_blob)) { ## innerj=1
        oldAPHLs <- .calc_APHLs_from_ZX(sXaug=sXaug, processed=processed, phi_est=phi_est, which=processed$p_v_obj, 
                                        lambda_est=lambda_est, dvdu=wranefblob$dvdu, u_h=u_h, mu=muetablob$mu)
      } 
      if ( ! GLMMbool) {
        # arguments for init_resp_z_corrections_new called in calc_zAug_not_LMM
        init_z_args <- c(constant_init_z_args,
                         list(w.ranef=wranefblob$w.ranef, u_h=u_h, v_h=v_h, dvdu=wranefblob$dvdu, 
                              sXaug=sXaug, w.resid=w.resid))
      } else init_z_args <- NULL
      calc_zAug_args <- c(constant_zAug_args,
                          list(muetablob=muetablob, dlogWran_dv_h=wranefblob$dlogWran_dv_h, 
                               sXaug=sXaug, 
                               w.ranef=wranefblob$w.ranef, 
                               w.resid=w.resid,
                               ############################# ZAL_scaling=ZAL_scaling,
                               init_z_args=init_z_args) )
      zInfo <- do.call(".calc_zAug_not_LMM",calc_zAug_args) 
      if (GLMMbool) zInfo$z2 <- NULL 
      etamo <- muetablob$sane_eta - off
      zInfo$z1_eta <- zInfo$z1-etamo
      z1_sscaled_eta <- zInfo$z1_sscaled - etamo # zAug[-seq_n_u_h]-etamo # z_1-sscaled-etamo
      if (GLMMbool) {
        zInfo$dlogfvdv <-  - v_h * wranefblob$w.ranef
      } else zInfo$dlogfvdv <- (zInfo$z2 - v_h) * wranefblob$w.ranef
      ## the gradient for -p_v (or -h), independent of the scaling
      zInfo$m_grad_obj <- c( ## drop() avoids c(Matrix..); attr(sXaug,"w.resid") presumably correct for truncated models too.
        m_grad_v <- drop(.crossprod(ZAL, attr(sXaug,"w.resid") * zInfo$z1_eta) + zInfo$dlogfvdv), # Z'W(z_1-eta)+ dlogfvdv
        drop(.crossprod(X.pv, attr(sXaug,"w.resid") * z1_sscaled_eta)) # X'W(z_1-sscaled-eta)
      )
      zInfo$gainratio_grad <- zInfo$m_grad_obj # needed for .do_damped_WLS_v_in_b_spprec()
      zInfo$scaled_grad <- zInfo$m_grad_obj # used as LM_z=zInfo["scaled_grad"] in .do_damped_WLS_v_in_b_spprec
      if (trace>1L) {
        maxs_grad <- c(max(abs(zInfo$m_grad_obj[seq_n_u_h])),max(abs(zInfo$m_grad_obj[-seq_n_u_h])))
        cat(stylefn(paste0("v_h iter=",innerj,", max(|grad|): v=",maxs_grad[1L],"beta=",maxs_grad[2L],";")))
      }
      constant_APHLs_args <- list(processed=processed, which="hlik", sXaug=sXaug, phi_est=phi_est, lambda_est=lambda_est)
      # the following block needs m_grad_v the new m_grad_v hence its position
      pot4improv <- get_from_MME(sXaug=sXaug, which="Mg_invH_g", B=m_grad_v)
      low_pot <- (pot4improv < pot_tol)
      damped_WLS_blob <- .do_damped_WLS_v_in_b_spprec(sXaug=sXaug, zInfo=zInfo, 
                                        old_Vscaled_beta=Vscaled_beta,
                                        oldAPHLs=oldAPHLs,
                                        APHLs_args = constant_APHLs_args,
                                        damping=.get_new_damping(dampings_env$v[["v_in_b"]],"v_in_b"),
                                        Trace=trace,
                                        ypos=ypos,off=off,
                                        GLMMbool=GLMMbool,etaFix=etaFix,
                                        constant_u_h_v_h_args=constant_u_h_v_h_args,
                                        updateW_ranefS_constant_arglist=updateW_ranefS_constant_arglist,
                                        wranefblob=wranefblob,seq_n_u_h=seq_n_u_h,
                                        update_sXaug_constant_arglist=update_sXaug_constant_arglist,
                                        ZAL_scaling=ZAL_scaling,
                                        processed=processed, 
                                        phi_est=phi_est, #H_global_scale=H_global_scale, 
                                        n_u_h=n_u_h, ZAL=ZAL,
                                        which_LevMar_step = "v",
                                        low_pot = structure(low_pot,pot_tol=pot_tol),
                                        stylefn=stylefn, # i.e., .spaMM.data$options$stylefns$vloop
                                        outer=FALSE
      ) 
      for (st in c("Vscaled_beta","wranefblob","v_h","u_h","muetablob",
                   "w.resid", 
                   "sXaug")) assign(st,damped_WLS_blob[[st]]) 
      if ( ! GLMMbool ) {
        #Xscal <- damped_WLS_blob$Xscal ## contains ZAL with new scaling, but weight_X is not applied since it is applied only locally in the mMatrix_method.
        ZAL_scaling <- damped_WLS_blob$ZAL_scaling # presumably the TAGged 1
      }
      
      #  At this point all return elements are updated as function of the latest Vscaled_beta.
      #  In particular We need muetablob and (if ! LMM) sXaug, hence a lot of stuff.
      #####
      beta_eta <- Vscaled_beta$beta_eta
      ##### assessment of convergence
      if ( low_pot ) { ## where we can see a gradient "small but measurable" for a pot4improv "really small".
        #if (damped_WLS_blob$APHLs[["hlik"]] < oldAPHLs[["hlik"]]) browser() 
        if (trace>1L) {cat(stylefn(" \u2713"))} # check mark
        breakcond <- "low_pot"
        break 
      } 
      # test on max(abs(m_grad_obj[seq_n_u_h]))<0.001 does not seem right ... typically FALSE, pushing the loop to maxit.
      relV <- v_h*sqrt(wranefblob$w.ranef)  ## convergence on v_h relative to sqrt(w.ranef)
      if (innerj>1L) {
        abs_d_relV <- abs(relV - old_relV) 
        not_moving <- ( mean(abs_d_relV) < d_relV_tol )
        if (not_moving) {
          breakcond <- "not_moving"
          break
        }
      }       
      # Updating for next iteration:
      oldAPHLs <- damped_WLS_blob$APHLs
      old_relV <- relV
      dampings_env$v[["v_in_b"]] <- damped_WLS_blob$damping
    } ################ E N D LOOP ##############
    if (innerj==maxit.mean) {
      if (trace) cat(crayon::red("!"))
      breakcond <- "maxit"
    }
    break_info$IRLS_breakcond <- breakcond
    break_info$maxs_grad <- c(max(abs(zInfo$m_grad_obj[seq_n_u_h])),max(abs(zInfo$m_grad_obj[-seq_n_u_h])))
    names(beta_eta) <- colnames(X.pv)
    RESU <- list(sXaug=sXaug, 
                 #fitted=fitted, 
                 #weight_X=weight_X, 
                 nobs=nobs, pforpv=pforpv, seq_n_u_h=seq_n_u_h, u_h=u_h, 
                 muetablob=muetablob, 
                 lambda_est=lambda_est,
                 phi_est=phi_est,
                 beta_eta=beta_eta, w.resid=w.resid, wranefblob=wranefblob, 
                 v_h=v_h, eta=muetablob$sane_eta, innerj=innerj,
                 break_info=break_info
    )
    return(RESU)
  }
