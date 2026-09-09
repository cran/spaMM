# corr algo has Xscal that can be p4m_ized without affecting $AUGI0_ZX  
# The original design matrices are provided by $AUGI0_ZX$X.pv  and ZAL argument
# .
# .
# .
.solve_v_h_IRLS <- # only for LevM && is_HL_1_1 # Despite the name, this expects augmented system with X block.
  function(ZAL, # the Z that gives eta and then the muetablob # and some more...
           # X.pv removed as argument bc apparently never different from processed$AUGI0_ZX$X.pv
           # This fn has p4m code, see .makeMatp4m() calls and 'is_p4m_H' 
           y, ## could be taken from processed ? 
           n_u_h, 
           lambda_est, off, maxit.mean, etaFix,
           wranefblob, processed,
           #
           # X.pv=processed$AUGI0_ZX$X.pv, # the X that gives eta and the weights: not perturbed in p4m
           ## supplement for ! LMM
           phi_est, 
           ## supplement for LevenbergM 
           w.resid=NULL, 
           ## supplement for LevenbergM
           beta_eta,
           ## "supplement for ! GLMM" (??): this is always used in this fn
           u_h, v_h,  
           #
           trace=FALSE,
           stylefn=.spaMM.data$options$stylefns$vloop,
           looseness,
           LevMarblob=NULL,
           dampings_env,
           damped_WLS_v_in_b_fn
  ) {
    X.pv <- processed$AUGI0_ZX$X.pv
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

    eta <- off + drop(ZAL %*% v_h + X.pv %*% beta_eta) # the p4m eta is still the one with the uncorrected matrices
    muetablob <- .muetafn(eta=eta,BinomialDen=processed$BinomialDen,processed=processed, phi_est=phi_est) 
    ## weight_X and Xscal varies within loop if ! LMM since at least the GLMweights in w.resid change
    w.resid <- .calc_w_resid(muetablob$GLMweights,phi_est, obsInfo=processed$how$obsInfo)
    # at this point w.resid is always the result of .calc_w_resid()
    # and when it is a list with info about mv model it has a complete vector $w_resid.
    H_w.resid <- .calc_H_w.resid(w.resid, muetablob=muetablob, processed=processed) # for LLF w.resid is not generally defined.
    H_global_scale <- .calc_H_global_scale(H_w.resid)
    ##### initial sXaug
    ZAL_scaling <- 1/sqrt(wranefblob$w.ranef*H_global_scale) ## Q^{-1/2}/s
    Xscal <- .make_Xscal(ZAL, ZAL_scaling = ZAL_scaling, processed=processed, as_matrix=.eval_as_mat_arg(processed))
    weight_X <- .calc_weight_X(Hobs_w.resid=H_w.resid, H_global_scale=H_global_scale, obsInfo=processed$how$obsInfo) ## sqrt(s^2 W.resid)  
    if (is_p4m_H <- ! is.null((multinom_info <- processed$multinom_info)[["mnsizes"]])) {
      p4mprobs <- .calc_p4mprobs(muetablob=muetablob, multinom_info)
      dcdv_p4m <- .makeMatp4m(mat=ZAL, multinom_info=multinom_info, processed=processed, 
                              p4mprobs=p4mprobs)
      dcdb_p4m <- .makeMatp4m(mat=X.pv, multinom_info=multinom_info, processed=processed, 
                              p4mprobs=p4mprobs)
      replaces_etamo <- drop(dcdv_p4m %*% v_h + dcdb_p4m %*% beta_eta)
      muetablob$dz1_p4m <- replaces_etamo -(eta-off)
      constant_zAug_args$ZAL <- dcdv_p4m # "doSeeMe" # dcdv_p4m seems logical    
      #             given this 'ZAL' is used for the y2_sscaled term, in factor with sscaled (\varsigma)
      #             so not from the term "in red". 
      Xscal_H_p4m <- .make_Xscal(dcdv_p4m, ZAL_scaling = ZAL_scaling, processed=processed, 
                                 as_matrix=.eval_as_mat_arg(processed), X=dcdb_p4m)
      sXaug <- do.call(processed$sXaug_method, # H_w.resid provided as attr(weight_X,"H_w.resid")!
                       list(Xaug=Xscal_H_p4m, weight_X=weight_X, w.ranef=wranefblob$w.ranef, H_global_scale=H_global_scale))
    } else sXaug <- do.call(processed$sXaug_method, # H_w.resid provided as attr(weight_X,"H_w.resid")!
                     list(Xaug=Xscal, weight_X=weight_X, w.ranef=wranefblob$w.ranef, H_global_scale=H_global_scale))
    if (trace) {
      tracechar <- ifelse(.BLOB(sXaug)$nonSPD,"!",".")
      cat(stylefn(tracechar)) # hmff blue (vloop) F I X M E
    }
    
    Vscaled_beta <- c(v_h/ZAL_scaling ,beta_eta)
    
    break_info <- list()
    which_i_llblock <- .which_i_llblock(Xscal, n_u_h) # preprocessing for faster updating of (sparse) Xscal when scaling changes
    
    ################ L O O P ##############
    
    for (innerj in 1:maxit.mean) {
      ##### get the lik of the current state
      if (is.null(damped_WLS_blob)) { ## innerj=1
        oldAPHLs <- .calc_APHLs_from_ZX(sXaug=sXaug, processed=processed, phi_est=phi_est, which=processed$p_v_obj, 
                                        lambda_est=lambda_est, dvdu=wranefblob$dvdu, u_h=u_h, muetablob=muetablob)
      }
      if ( ! GLMMbool) {
        # arguments for init_resp_z_corrections_new called in calc_zAug_not_LMM
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
      ## the gradient for -p_v (or -h), independent of the scaling
      if (zInfo$z1_is4p4m) { 
        dcdmu <- zInfo$z1-replaces_etamo 
      } else {
        etamo <- muetablob$sane_eta-off
        dcdmu <- zInfo$z1-etamo
      }
      if (is_p4m_H) {
        m_grad_obj <- .calc_m_grad_obj(zInfo, dcdmu=dcdmu, GLMMbool=GLMMbool, v_h=v_h, 
                                       wranefblob=wranefblob, H_w.resid=.BLOB(sXaug)$H_w.resid, 
                                       dLinkPred_dv=dcdv_p4m, dLinkPred_db=dcdb_p4m)
      } else m_grad_obj <- .calc_m_grad_obj(zInfo, dcdmu=dcdmu, GLMMbool=GLMMbool, v_h=v_h, 
                                            wranefblob=wranefblob, H_w.resid=.BLOB(sXaug)$H_w.resid, 
                                            dLinkPred_dv=ZAL, dLinkPred_db=X.pv)
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
      damped_WLS_blob <- 
        damped_WLS_v_in_b_fn( # performs a loop of damping_to_solve() calls until some improvement of objfn.
          sXaug=sXaug, zInfo=zInfo, ZAL=ZAL,
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
          phi_est=phi_est, H_global_scale=H_global_scale, n_u_h=n_u_h, 
          which_i_llblock=which_i_llblock,
          which_LevMar_step = "v", # which implies that sXaug is not updated in each iteration of the called loop, 
                                   # but it may be updated at its end.
          low_pot = structure(low_pot,pot_tol=pot_tol),
          stylefn=stylefn, # i.e., .spaMM.data$options$stylefns$vloop
          outer=FALSE) 
      list2env(damped_WLS_blob[c("w.resid", ## !important! cf test-adjacency-corrMatrix.R
                                 "Vscaled_beta","wranefblob","v_h","u_h","muetablob", "weight_X", "sXaug",
                                 "dcdv_p4m")], 
               envir = environment()) 
      # for (st in c("Vscaled_beta","wranefblob","v_h","u_h","muetablob",
      #              "w.resid", ## !important! cf test-adjacency-corrMatrix.R
      #              "weight_X", 
      #              "sXaug")) assign(st,damped_WLS_blob[[st]]) 
      if (is_p4m_H) { 
        constant_zAug_args$ZAL <- dcdv_p4m # "doSeeMe" # dcdv_p4m seems logical    
        #   #             given this 'ZAL' is used for the y2_sscaled term, in factor with sscaled (\varsigma)
        #   #             so not from the term "in red". 
      } else if ( ! GLMMbool ) {
        Xscal <- damped_WLS_blob$Xscal ## contains ZAL with new scaling, but weight_X is not applied since it is applied only locally in the sXaug_method
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
      if (trace) cat(cli::col_red("!")) ## _F I X M E_ I could have tried to assess the speed of convergence in order to decide whether to exit the loop or not...
      breakcond <- "maxit"
    }
    break_info$IRLS_breakcond <- breakcond
    if (pforpv) { 
      break_info$maxs_grad <- c(max(abs(m_grad_obj[seq_n_u_h])),max(abs(m_grad_obj[-seq_n_u_h])))
    } else break_info$maxs_grad <- c(max(abs(m_grad_obj[seq_n_u_h])), 0)
    names(beta_eta) <- colnames(X.pv)
    if (! is.null(damped_WLS_blob)) {
      fitted <- damped_WLS_blob$fitted
      weight_X <- damped_WLS_blob$weight_X ## F I X M E it seems better to store  weight_X  as attr(sXaug,...) and no weight_X elsewhere in output 
    } 
    RESU <- list(# sXaug=sXaug, 
                 # fitted=fitted, 
                 # weight_X=weight_X, 
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
