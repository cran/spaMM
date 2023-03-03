.solve_v_h_IRLS_spprec <- # only for LevM && is_HL_1_1
  function(X.pv, 
           ZAL, y, ## could be taken fom processed ? 
           n_u_h, 
           #H_global_scale, 
           lambda_est, off, maxit.mean, etaFix,
           wranefblob, processed,
           ## supplement for ! LMM
           phi_est, 
           ## supplement for for LevenbergM 
           w.resid=NULL, 
           ## supplement for LevenbergM
           beta_eta,
           ## supplement for ! GLMM (??)
           u_h, v_h,  
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
    LMMbool <-  attr(processed[["models"]],"LMMbool")
    GLMMbool <- attr(processed[["models"]],"GLMMbool") 
    not_moving <- FALSE
    old_relV <- NULL
    damped_WLS_blob <- NULL
    # looseness controlled by .wrap_do_damped_WLS_outer(..., looseness, ...) which may be overwritten when "strict_v|b"
    pot_tol <- processed$spaMM_tol$v_pot_tol * looseness
    d_relV_tol <- processed$spaMM_tol$d_relV_tol * looseness
    constant_zAug_args <- list(n_u_h=n_u_h, nobs=nobs, pforpv=pforpv, y=y, off=off, ZAL=ZAL, processed=processed)

    ##### initial sXaug
    ZAL_scaling <- 1  ## TAG: scaling for spprec
    eta  <- off + drop(processed$AUGI0_ZX$X.pv %*% beta_eta) + drop(ZAL %id*% v_h) 
    muetablob <- .muetafn(eta=eta,BinomialDen=processed$BinomialDen,processed=processed, phi_est=phi_est) 
    ## weight_X and Xscal varies within loop if ! LMM since at least the GLMweights in w.resid change
    w.resid <- .calc_w_resid(muetablob$GLMweights,phi_est, obsInfo=processed$how$obsInfo)
    ## needs adjMatrix and corrPars to define Qmat
    update_sXaug_constant_arglist <- list(AUGI0_ZX=processed$AUGI0_ZX, corrPars=corrPars, cum_n_u_h=processed$cum_n_u_h) 
    sXaug_arglist <- c(update_sXaug_constant_arglist,
                       list(w.ranef=wranefblob$w.ranef,
                            H_w.resid=.calc_H_w.resid(w.resid, muetablob=muetablob, processed=processed)))
    # .HLfit_body_augZXy has called .init_AUGI0_ZX_envir_spprec_info(processed,LMatrices)...
    sXaug <- do.call(processed$spprec_method, # ie, def_AUGI0_ZX_spprec
                     sXaug_arglist)
    if (trace) {
      tracechar <- ifelse(.BLOB(sXaug)$nonSPD,"!",".")
      cat(stylefn(tracechar)) # hmff blue (vloop) F I X M E
    }
    
    Vscaled_beta <- list(v_h=v_h/ZAL_scaling,beta_eta=beta_eta)
    
    break_info <- list()
    ################ L O O P ##############
    for (innerj in 1:maxit.mean) {
      ##### get the lik of the current state
      if (is.null(damped_WLS_blob)) { ## innerj=1
        oldAPHLs <- .calc_APHLs_from_ZX(sXaug=sXaug, processed=processed, phi_est=phi_est, which=processed$p_v_obj, 
                                        lambda_est=lambda_est, dvdu=wranefblob$dvdu, u_h=u_h, muetablob=muetablob)
      } 
      if ( ! GLMMbool) {
        # # arguments for init_resp_z_corrections_new called in calc_zAug_not_LMM
        # init_z_args <- c(constant_init_z_args,
        #                  list(w.ranef=wranefblob$w.ranef, u_h=u_h, v_h=v_h, dvdu=wranefblob$dvdu, 
        #                       sXaug=sXaug))  # H_w.resid provided by sXaug!
        # z2 <- do.call(".init_resp_z_corrections_new",init_z_args)$z20
        z2 <- .calc_z2(lcrandfamfam=lcrandfamfam, psi_M=processed$psi_M, cum_n_u_h=processed$cum_n_u_h, rand.families=processed$rand.families, 
                       u_h=u_h, lambda_est=lambda_est, v_h=v_h, dvdu=wranefblob$dvdu)
      } else z2 <- rep(0,n_u_h)
      calc_zAug_args <- c(constant_zAug_args,
                          list(muetablob=muetablob, dlogWran_dv_h=wranefblob$dlogWran_dv_h, 
                               sXaug=sXaug,  # H_w.resid provided by sXaug!
                               w.ranef=wranefblob$w.ranef, 
                               w.resid=w.resid,
                               z2=z2) )
      zInfo <- do.call(".calc_zAug_not_LMM",calc_zAug_args) 
      if (GLMMbool) zInfo$z2 <- NULL 
      zInfo$m_grad_obj <- .calc_m_grad_obj(zInfo, GLMMbool=GLMMbool, v_h=v_h, 
                                     wranefblob=wranefblob, H_w.resid=.BLOB(sXaug)$H_w.resid, ZAL=ZAL, X.pv=X.pv, etamo=muetablob$sane_eta - off)

      zInfo$gainratio_grad <- zInfo$m_grad_obj # needed for .do_damped_WLS_v_in_b_spprec()
      zInfo$scaled_grad <- zInfo$m_grad_obj # used as LM_z=zInfo["scaled_grad"] in .do_damped_WLS_v_in_b_spprec
      if (trace>1L) {
        if (pforpv) { 
          maxs_grad <- c(max(abs(zInfo$m_grad_obj[seq_n_u_h])),max(abs(zInfo$m_grad_obj[-seq_n_u_h])))
        } else maxs_grad <- c(max(abs(zInfo$m_grad_obj[seq_n_u_h])),0) # outer beta
        cat(stylefn(paste0("v_h iter=",innerj,", max(|grad|): v=",maxs_grad[1L],"beta=",maxs_grad[2L],";")))
      }
      constant_APHLs_args <- list(processed=processed, which="hlik", sXaug=sXaug, phi_est=phi_est, lambda_est=lambda_est)
      # the following block needs m_grad_v the new m_grad_v hence its position
      m_grad_v <- zInfo$m_grad_obj[seq_n_u_h]
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
                                        lambda_est=lambda_est,
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
      list2env(damped_WLS_blob[c("w.resid", ## !important! 
                                 "Vscaled_beta","wranefblob","v_h","u_h","muetablob", "sXaug")], envir = environment()) 
      # for (st in c("Vscaled_beta","wranefblob","v_h","u_h","muetablob",
      #              "w.resid", 
      #              "sXaug")) assign(st,damped_WLS_blob[[st]]) 
      if ( ! GLMMbool ) {
        #Xscal <- damped_WLS_blob$Xscal ## contains ZAL with new scaling, but weight_X is not applied since it is applied only locally in the corr_method
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
    if (pforpv) { 
      break_info$maxs_grad <- c(max(abs(zInfo$m_grad_obj[seq_n_u_h])),max(abs(zInfo$m_grad_obj[-seq_n_u_h])))
    } else break_info$maxs_grad <- c(max(abs(zInfo$m_grad_obj[seq_n_u_h])),0) # outer beta
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
