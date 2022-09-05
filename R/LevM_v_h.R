.solve_v_h_IRLS <- # only for LevM && is_HL_1_1
  function(X.pv, 
           ZAL, y, ## could be taken fom processed ? 
           n_u_h, 
           H_global_scale, 
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
    if ( ! GLMMbool) {
      constant_init_z_args <- c(list(lcrandfamfam=lcrandfamfam, nobs=nobs, lambda_est=lambda_est, ZAL=ZAL),  
                                # fit_as_ZX args specific for ! GLMM:
                                for_init_z_args,
                                #
                                mget(c("cum_n_u_h","rand.families"),envir=processed))
    } 
    
    ##### initial sXaug
    ZAL_scaling <- 1/sqrt(wranefblob$w.ranef*H_global_scale) ## Q^{-1/2}/s
    Xscal <- .make_Xscal(ZAL, ZAL_scaling = ZAL_scaling, processed=processed, as_matrix=.eval_as_mat_arg(processed))
    if (inherits(Xscal,"Matrix")) { # same type as ZAL
      #@p[c] must contain the index _in @x_ of the first nonzero element of column c, x[p[c]] in col c and row i[p[c]])  
      elmts_affected_cols <- seq_len(Xscal@p[n_u_h+1L]) ## corresponds to cols seq_n_u_h
      which_i_affected_rows <- which(Xscal@i[elmts_affected_cols]>(n_u_h-1L))    
    } else {
      which_i_affected_rows <- NULL
    }
    if (is.null(muetablob)) { ## NULL input eta allows NULL input muetablob
      eta  <- off + (Xscal %*% c(v_h/ZAL_scaling ,beta_eta))[ypos]
      muetablob <- .muetafn(eta=eta,BinomialDen=processed$BinomialDen,processed=processed) 
    }
    ## weight_X and Xscal varies within loop if ! LMM since at least the GLMweights in w.resid change
    if ( is.null(w.resid) ) w.resid <- .calc_w_resid(muetablob$GLMweights,phi_est)
    # at this point w.resid is always the result of .calc_w_resid()
    # and when it is a list with info about mv model it has a complete vector $w_resid.
    H_w.resid <- .calc_H_w.resid(w.resid, muetablob=muetablob, processed=processed) # for LLF w.resid is not generally defined.
    weight_X <- .calc_weight_X(Hobs_w.resid=H_w.resid, H_global_scale=H_global_scale, obsInfo=processed$how$obsInfo) ## sqrt(s^2 W.resid)  
    sXaug <- do.call(processed$mMatrix_method,
                     list(Xaug=Xscal, weight_X=weight_X, w.ranef=wranefblob$w.ranef, H_global_scale=H_global_scale))
    if (trace) {
      tracechar <- ifelse(identical(attr(sXaug,"BLOB")$nonSPD,TRUE),"!",".")
      cat(stylefn(tracechar)) # hmff blue (vloop) F I X M E
    }
    
    Vscaled_beta <- c(v_h/ZAL_scaling ,beta_eta)
    
    break_info <- list()
    ################ L O O P ##############
    for (innerj in 1:maxit.mean) {
      ##### get the lik of the current state
      if (is.null(damped_WLS_blob)) { ## innerj=1
        oldAPHLs <- .calc_APHLs_from_ZX(sXaug=sXaug, processed=processed, phi_est=phi_est, which=processed$p_v_obj, 
                                        lambda_est=lambda_est, dvdu=wranefblob$dvdu, u_h=u_h, muetablob=muetablob)
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
      wzAug <- c(zInfo$y2_sscaled/ZAL_scaling, (zInfo$z1_sscaled)*weight_X) 
      etamo <- muetablob$sane_eta - off
      ## the gradient for -p_v (or -h), independent of the scaling
      m_grad_obj <- .calc_m_grad_obj(zInfo, z1_eta=zInfo$z1-etamo, z1_sscaled_eta=zInfo$z1_sscaled - etamo, GLMMbool, v_h, wranefblob, 
                                     H_w.resid=.BLOB(sXaug)$H_w.resid, ZAL, X.pv)
      ## amounts to 
      # m_grad_obj <- c( ## drop() avoids c(Matrix..) 
      #   m_grad_v <- drop(.crossprod(ZAL, WLS_mat_weights * zInfo$z1_eta) + dlogfvdv), # Z'W(z_1-eta)+ dlogfvdv 
      #   drop(.crossprod(X.pv, WLS_mat_weights * z1_sscaled_eta)) # X'W(z_1-sscaled-eta)
      # )
      ## where the zInfo terms must depend on WLS_mat_weights too! 
      if (trace>1L) {
        if (pforpv) { 
          maxs_grad <- c(max(abs(m_grad_obj[seq_n_u_h])),max(abs(m_grad_obj[-seq_n_u_h])))
        } else maxs_grad <- c(max(abs(m_grad_obj[seq_n_u_h])), 0) # outer beta
        cat(stylefn(paste0("v_h iter=",innerj,", max(|grad|): v=",maxs_grad[1L],"beta=",maxs_grad[2L],";")))
      }
      zInfo$gainratio_grad <- m_grad_obj ## before rescaling
      # gradient for scaled system from gradient of objective
      scaled_grad <- H_global_scale * m_grad_obj
      scaled_grad[seq_n_u_h] <- scaled_grad[seq_n_u_h] * ZAL_scaling 
      zInfo$scaled_grad <- scaled_grad
      constant_APHLs_args <- list(processed=processed, which="hlik", sXaug=sXaug, phi_est=phi_est, lambda_est=lambda_est)
      # the following block needs m_grad_v the new m_grad_v hence its position
      m_grad_v <- m_grad_obj[seq_n_u_h]
      pot4improv <- get_from_MME(sXaug=sXaug, which="Mg_invH_g", B=m_grad_v)
      low_pot <- (pot4improv < pot_tol)
      damped_WLS_blob <- .do_damped_WLS_v_in_b(sXaug=sXaug, zInfo=zInfo, 
                                        old_Vscaled_beta=Vscaled_beta,
                                        oldAPHLs=oldAPHLs,
                                        APHLs_args = constant_APHLs_args,
                                        damping=.get_new_damping(dampings_env$v[["v_in_b"]],"v_in_b"),
                                        Trace=trace,
                                        ypos=ypos,off=off,
                                        GLMMbool=GLMMbool,etaFix=etaFix,
                                        lambda_est=lambda_est,
                                        wranefblob=wranefblob,seq_n_u_h=seq_n_u_h,ZAL_scaling=ZAL_scaling,
                                        processed=processed, Xscal=Xscal,
                                        phi_est=phi_est, H_global_scale=H_global_scale, n_u_h=n_u_h, ZAL=ZAL,
                                        which_i_affected_rows=which_i_affected_rows,
                                        which_LevMar_step = "v",
                                        low_pot = structure(low_pot,pot_tol=pot_tol),
                                        stylefn=stylefn, # i.e., .spaMM.data$options$stylefns$vloop
                                        outer=FALSE
      ) 
      list2env(damped_WLS_blob[c("w.resid", ## !important! cf test-adjacency-corrMatrix.R
                                 "Vscaled_beta","wranefblob","v_h","u_h","muetablob", "weight_X", "sXaug")], envir = environment()) 
      # for (st in c("Vscaled_beta","wranefblob","v_h","u_h","muetablob",
      #              "w.resid", ## !important! cf test-adjacency-corrMatrix.R
      #              "weight_X", 
      #              "sXaug")) assign(st,damped_WLS_blob[[st]]) 
      if ( ! GLMMbool ) {
        Xscal <- damped_WLS_blob$Xscal ## contains ZAL with new scaling, but weight_X is not applied since it is applied only locally in the mMatrix_method.
        ZAL_scaling <- damped_WLS_blob$ZAL_scaling
      }
      
      #  At this point all return elements are updated as function of the latest Vscaled_beta.
      #  In particular We need muetablob and (if ! LMM) sXaug, hence a lot of stuff.
      #####
      beta_eta <- Vscaled_beta[n_u_h+seq_len(pforpv)]
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
      if (trace) cat(crayon::red("!")) ## _F I X M E_ I could have tried to assess the speed of convergence in order to decide whether to exit the loop or not...
      breakcond <- "maxit"
    }
    break_info$IRLS_breakcond <- breakcond
    if (pforpv) { # outer beta: several change in v3.11.5 without parallel changes in .solve_IRLS_as_spprec
      break_info$maxs_grad <- c(max(abs(m_grad_obj[seq_n_u_h])),max(abs(m_grad_obj[-seq_n_u_h])))
    } else break_info$maxs_grad <- c(max(abs(m_grad_obj[seq_n_u_h])), 0)
    names(beta_eta) <- colnames(X.pv)
    if (! is.null(damped_WLS_blob)) {
      fitted <- damped_WLS_blob$fitted
      weight_X <- damped_WLS_blob$weight_X ## F I X M E it seems better to store  weight_X  as attr(sXaug,...) and no weight_X elsewhere in output 
    } 
    RESU <- list(sXaug=sXaug, 
                 fitted=fitted, 
                 weight_X=weight_X, 
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
