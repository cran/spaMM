.do_damped_WLS <- function(sXaug, # constant within the fn, which instead evaluates newsXaug and RESU$sXaug.
                           zInfo, # fixed descriptors of system to be solved
           old_Vscaled_beta,
           oldAPHLs,
           APHLs_args,
           damping,  # _F I X M E_ looks like it could be parallelized 
           dampingfactor=2, ## no need to change the input value
           ypos, off, GLMMbool, etaFix, 
           lambda_est,
           wranefblob, seq_n_u_h, 
           ZAL_scaling,
           processed, 
           Trace, 
           phi_est, H_global_scale, n_u_h, 
           ZAL,
           which_LevMar_step, # fixed within this function
           which_i_llblock,
           Xscal, ## locally fixed, "resident"
           corr_method=processed$corr_method,
           # promise rather than argument:
           low_pot=NULL,
           v_infer_args=NULL, # not null for beta optimization with v_in_b optimization i.e. in *some* .do_damped_WLS_outer() call
           stylefn, # in-loop stylefn for damped_WLS
           stylefn_v_out= .spaMM.data$options$stylefns$v_out_last, 
           stylefn_v_in= .spaMM.data$options$stylefns$v_in_last, ##
           outer,
           IRLS_fn=get(".solve_v_h_IRLS", asNamespace("spaMM"), inherits=FALSE) 
) {
  if (outer) {
    trace <- max(0L,Trace-1L)
    stylefn_v <- stylefn_v_out
  } else {
    trace <- max(0L,Trace-2L) # => TRACE=3 will give more details on steps of v_h_IRLS thant TRACE=2
    stylefn_v <- stylefn_v_in
  }
  if (processed$p_v_obj=="p_v" && which_LevMar_step!="v") { 
    objname <- "p_v" 
  } else { 
    objname <- APHLs_args$which <- "hlik"
  }
  oldlik <- oldAPHLs[[objname]]
  newXscal <- Xscal ## template
  initdamping <- damping
  gainratio_grad <- zInfo$gainratio_grad
  # grad wrt scaled v = d f / d (v/ZAL_scaling) = ZAL_scaling * d f / d v
  restarted_at_e_7 <- FALSE
  first_it <- TRUE
  prev_gainratio <- -Inf
  if (Trace && ! is.null(v_infer_args)) {
    cat(stylefn("[")) # cat(which_LevMar_step) #=> a substep of V_IN_B: "strict_v|b" or "b_&_v_in_b"
  } 
  GLGLLM_const_w <- attr(processed$models,"GLGLLM_const_w")
  while ( TRUE ) { ## loop on damping; each iteration produce blue + ul-greens + yellow
    if (processed$HL[1L]==1L) { ## ML fit 
      Vscaled_beta <- old_Vscaled_beta
      ## maximize p_v wrt beta only
      if ( which_LevMar_step=="v_b") { ## note tests on summand too !!!!!!!
        LevMarblob <- get_from_MME(sXaug=sXaug, which="LevMar_step", LMrhs=zInfo$scaled_grad, damping=damping)
        Vscaled_beta <- Vscaled_beta + LevMarblob$dVscaled_beta 
      } else if ( which_LevMar_step=="v") { ## v_h estimation given beta 
        LevMarblob <- get_from_MME(sXaug=sXaug, which="LevMar_step_v_h", LMrhs=zInfo$scaled_grad[seq_n_u_h], damping=damping)
        Vscaled_beta[seq_n_u_h] <- Vscaled_beta[seq_n_u_h] + LevMarblob$dVscaled 
      } else if ( which_LevMar_step=="strict_v|b") {
        # =: the case where we only fit v_h for the input beta_eta:
        #    only v_h in Vscaled_beta will be changed, by v_infer_args step below, 
        LevMarblob <- v_infer_args$LevMarblob ## LevMarblob$dVscaled_beta will be overwritten below
        #                                        is.null(LevMarblob) occurs in initial damping=Inf call.
      } else if ( which_LevMar_step=="b_&_v_in_b") { 
        LevMarblob <- get_from_MME(sXaug=sXaug, which="LevMar_step", LMrhs=zInfo$scaled_grad, damping=damping) # template
        dbeta_LevMarblob <- get_from_MME(sXaug=sXaug, which="LevMar_step_beta", LMrhs=zInfo$scaled_grad[-seq_n_u_h], damping=damping) # template
        Vscaled_beta[-seq_n_u_h] <- Vscaled_beta[-seq_n_u_h] + dbeta_LevMarblob$dbeta
        # v_h in Vscaled_beta will be changed by v_infer_args step below.
        v_infer_args$LevMarblob <- get_from_MME(sXaug=sXaug, which="LevMar_step_v_h", LMrhs=zInfo$scaled_grad[seq_n_u_h], damping=damping)
      } else if ( which_LevMar_step=="b_from_v_b") { 
        LevMarblob <- get_from_MME(sXaug=sXaug, which="LevMar_step", LMrhs=zInfo$scaled_grad, damping=damping)
        Vscaled_beta[-seq_n_u_h] <- Vscaled_beta[-seq_n_u_h] + LevMarblob$dVscaled_beta[-seq_n_u_h]
      } else if ( which_LevMar_step=="b") { ## currently not used ## If called when pforpv=0, get_from_MME -> .damping_to_solve() crashes badly.
        LevMarblob <- get_from_MME(sXaug=sXaug, which="LevMar_step_beta", LMrhs=zInfo$scaled_grad[-seq_n_u_h], damping=damping)
        Vscaled_beta[-seq_n_u_h] <- Vscaled_beta[-seq_n_u_h] + LevMarblob$dbeta 
      } 
      if ( which_LevMar_step!="v" &&  ! is.null(v_infer_args)) {
        
        if (Trace) cat(stylefn_v("["))
        v_h_blob <- .wrap_wrap_v_h_IRLS(IRLS_fn=IRLS_fn, v_h=Vscaled_beta[seq_n_u_h] * ZAL_scaling, 
                                        beta_eta=Vscaled_beta[-seq_n_u_h], seq_n_u_h, GLMMbool, wranefblob, processed, 
                                        lambda_est, v_infer_args, Trace) # each underline green is a damping _loop_ not a damping step
        if (Trace) cat(stylefn_v("]"))
        if (trace) .cat_break_info( v_h_blob, stylefn_v, stylefn)  ## prints, at the level of the outer damped_WLS, the results of the v_h IRLS
        
        Vscaled_beta[seq_n_u_h] <- v_h_blob$v_h/ZAL_scaling 
        LevMarblob$dVscaled_beta <- Vscaled_beta - old_Vscaled_beta
      }
    } else { ## joint hlik maximization
      LevMarblob <- get_from_MME(sXaug=sXaug, which="LevMar_step", LMrhs=zInfo$scaled_grad, damping=damping)
      Vscaled_beta <- old_Vscaled_beta + LevMarblob$dVscaled_beta 
    }
    fitted <- drop(Xscal %*% Vscaled_beta) ## length nobs+nr ! 
    eta <- fitted[ypos] + off
    newmuetablob <- .muetafn(eta=eta,BinomialDen=processed$BinomialDen,processed=processed, phi_est=phi_est) 
    fitted[ypos] <- newmuetablob$sane_eta
    
    neww.resid <- .calc_w_resid(newmuetablob$GLMweights,phi_est, obsInfo=processed$how$obsInfo)
    # at this point w.resid is always the result of .calc_w_resid()
    # and when it is a list with info about mv model it has a complete vector $w_resid.
    newH_w.resid <- .calc_H_w.resid(neww.resid, muetablob=newmuetablob, processed=processed) # for LLF w.resid is not generally defined.
    newweight_X <- .calc_weight_X(Hobs_w.resid=newH_w.resid, H_global_scale=H_global_scale, obsInfo=processed$how$obsInfo) ## sqrt(s^2 W.resid)  

    if (is.null(etaFix$v_h)) { 
      v_h <- Vscaled_beta[seq_n_u_h] * ZAL_scaling ## use original scaling!
      if (GLMMbool) {
        u_h <- v_h 
        newwranefblob <- wranefblob ## keep input wranefblob since GLMM and lambda_est not changed
      } else {
        u_h <- processed$u_h_v_h_from_v_h(v_h)
        if ( ! is.null(maybe <- attr(u_h,"v_h"))) v_h <- maybe
        ## update functions u_h,v_h
        newwranefblob <- wranefblob <- processed$updateW_ranefS(u_h=u_h,v_h=v_h, lambda=lambda_est)
      } 
    } else newwranefblob <- wranefblob
    mMatrix_arglist <- list(weight_X=newweight_X, w.ranef=newwranefblob$w.ranef, H_global_scale=H_global_scale)
    if ( ! GLMMbool ) {
      # newZAL_scaling necessary to get the correct logdet_sqrt_d2hdv2 for newsXaug
      newZAL_scaling <- 1/sqrt(newwranefblob$w.ranef*H_global_scale) ## Q^{-1/2}/s
      ## used only to compute a likelihood, not to update a system to be solved.
      mMatrix_arglist$Xaug <- newXscal <- .calc_Xscal_newscaled(newXscal, newZAL_scaling, ZAL, which_i_llblock, 
                                                    n_u_h, seq_n_u_h, processed)
    } else mMatrix_arglist$Xaug <- Xscal ##not distinct from the 'resident' Xscal
    ####
    APHLs_args$dvdu <- newwranefblob$dvdu
    APHLs_args$u_h <- u_h 
    APHLs_args$muetablob <- newmuetablob
    #
    if (processed$p_v_obj=="p_v" && which_LevMar_step!="v") { 
      if (GLGLLM_const_w) {
        newsXaug <- NULL
        APHLs_args$sXaug <- sXaug
      } else {
        newsXaug <- do.call(corr_method, mMatrix_arglist)
        tracechar <- ifelse(.BLOB(newsXaug)$nonSPD,"!",".")
        if (Trace) cat(stylefn(tracechar)) # yellow in V_IN_B case
        APHLs_args$sXaug <- newsXaug
      }
    } else { ## mMatrix_arglist will still be used after the loop !!!!!!!!!!!!!!!!!!!!
      APHLs_args$sXaug <- newsXaug <- NULL
    }
    newAPHLs <- do.call(".calc_APHLs_from_ZX", APHLs_args)
    #print(c(unlist(newAPHLs)))
    newlik <- unlist(newAPHLs[objname]) # keep name
    #
    if (damping==0L) {
      breakcond <- "damping=0"
      break
    }
    if ( which_LevMar_step=="strict_v|b") {
      breakcond <- "v|b_no_loop"
      attr(breakcond,"v_pot4improv") <- .pot4improv("v", sXaug, gainratio_grad=zInfo$gainratio_grad, seq_n_u_h)
      break
    } #  =: single call to .calc_APHLs_from_ZX to only fit v_h for the input beta_eta.
    if (first_it) { # test run many times, may be true only the first time
      pot4improv <- .pot4improv(which_LevMar_step, sXaug, gainratio_grad=zInfo$gainratio_grad, seq_n_u_h)
      loc_pot_tol <- .loc_pot_tol(which_LevMar_step, processed$spaMM_tol)
      if (is.null(low_pot)) low_pot <- (pot4improv < loc_pot_tol) 
    } # else ... would be highly suspect here; at this point, loc_pot_tol (necess to assess "stuck_obj") must be available
    if (low_pot) { # keeping the low_pot condition may be important for the "605" tests. We may always suppress it by the spaMM.options.
      very_low_pot <- (pot4improv < loc_pot_tol/10)
      breakcond <- "low_pot"
      attr(breakcond, "pot4improv") <- pot4improv
      attr(breakcond, "very_low_pot") <- very_low_pot
      if (processed$HL[1L]==1L && which_LevMar_step=="v_b") {
        v_pot4improv <- .pot4improv("v", sXaug, gainratio_grad=zInfo$gainratio_grad, seq_n_u_h)
        attr(breakcond, "no_overfit") <- ((v_pot4improv < processed$spaMM_tol$v_pot_tol))
      } 
      break
    } 
    ## ELSE
    gainratio <- (newlik!=-Inf) ## -Inf occurred in binary probit with extreme eta... and playing with hard COMPoisson cases...
    if (gainratio) { # non-zero (i.e., including negative values)
      summand <- .calc_summand_gainratio(processed, which_LevMar_step, LevMarblob, seq_n_u_h, ZAL_scaling, gainratio_grad)
      ## The two terms of the summand should be positive. In part. conv_dbetaV*rhs should be positive. 
      ## However, numerical error may lead to <0 or even -Inf
      ## Further, if there are both -Inf and +Inf elements the sum is NaN.
      summand[summand<0] <- 0
      denomGainratio <- sum(summand)
      dlogL <- newlik-oldlik
      conv_logL <- abs(dlogL)/(1+abs(newlik))
      gainratio <- 2*dlogL/denomGainratio ## cf computation in MadsenNT04 below 3.14, but generalized for D' != I ; numerator is reversed for maximization
      # but if gradient is practically zero and damping  ~0 we may not wish to compare ~0 to ~0...
      #  which is why we break, rather than stop, on (damping>1e100). MadsenNT04 have other stopping crits 
    } else { ## 2017/10/16 a patch to prevent a stop in case gainratio=0, but covers up dubious computations (FIXME)
      newlik <- -.Machine$double.xmax
      dlogL <- newlik-oldlik
      conv_logL <- abs(dlogL)/(1+abs(newlik))
      denomGainratio <- Inf # bc denomGainratio can be tested below 
    }
    if (trace) {
      cat(stylefn(paste(" dampingfactor=",dampingfactor,#"innerj=",innerj,
                           "damping=",damping,"gainratio=",gainratio,"oldlik=",oldlik,"newlik=",newlik))) # oldlik is from before the damping loop
    }
    if (is.nan(gainratio)) {
      # if the initial logL is the solution logL, then damping diverges 
      # it is then possible that some element of dVscaled_beta =0 and some of dampDpD =Inf
      # then the summand has some NaN
      # At the same time not all elements of dVscaled_beta need be 0 (eg small changes in eta for mu=0 or 1 in binomial models)
      # so testing dVscaled_beta is not sufficient to stop the algo
      # (LevenbergM is quite slow in such cases)
      breakcond <- "NaN_gain"
      break
    }
    if (objname == "p_v") { # then the starting objective may be 'too high' and we need to handle that
      div_gainratio <- (gainratio-prev_gainratio)*damping/dampingfactor
      if (div_gainratio < -max(LevMarblob$dampDpD)/(1e05*damping)) { ## 
        # tests: default tests + test_Materndifficult + test-nloptr
        breakcond <- "div_gain" # used to switch from "V_IN_B" to "strict_v|b"
        break 
      } 
    }
    if (gainratio > 0) { # gainratio may be negative if initial ranefs better optimize logL than the correct solution does.
      ## cf Madsen-Nielsen-Tingleff again, and as in levmar library by Lourakis
      damping <- damping * max(1/3,1-(2*gainratio-1)^3) # gainratio->0 factor-> 2; gainratio->1 factor->1/3
      breakcond <- "OK_gain"
      break 
    } 

    if (
      conditionNS_for_restart <-  (
        initdamping>1e-7 &&
        ( ! restarted_at_e_7 ) &&
        denomGainratio<loc_pot_tol/10
      )
    ) { # then restart 
      damping <- 1e-7
      dampingfactor <- 2
      restarted_at_e_7 <- TRUE
      prev_gainratio <- -Inf
      if (trace) cat("-")
      # and continue 
    } else if ( restarted_at_e_7 && ## has gone through 1e-7, # handles the cases where we started with too high damping 
                (( odd_condition <- (dampingfactor>4 && ## ie at least 2 iteration of the while() => prev_conv_logL is available
                                     conv_logL <1e-8 && abs(prev_conv_logL) <1e-8)) || # possible reason for odd_condition is truely stuck because of "COMPoisson" reasons?
                 damping>initdamping)
    ) { # restarted && damping>init must be a sufficient condition for break
      breakcond <- "stuck_obj"
      break   ##   cases were we do not expect any significant improvement
    } else { ## other UNsuccessful step
      prev_gainratio <- gainratio
      prev_conv_logL <- conv_logL
      damping <- damping*dampingfactor
      dampingfactor <- dampingfactor*2
      if (damping>1e10) { # endpoint for large negative gainratio (i.e. v overfit as starting point) => ___F I X M E___ do something clever in that case ?
        breakcond <- "div_damp"
        break 
      }
    }
    first_it <- FALSE # : skipped if break in first iteration
  } ################# end while(TRUE)
  if (Trace && ! is.null(v_infer_args)) {
    cat(stylefn("]"))
    if (trace) {cat(stylefn(damping))}
  }
  if (trace) cat(breakcond)
  if (is.null(newsXaug)) { ## which means that hlik is the local objective or that (GLGLLM_const_w).
    # For HL11, p_v will be used as oldAPHLs in the next call to .do_damped_WLS_outer() in an alternating algo;
    #   and sXaug may be needed to compute sscaled in .solve_v_h_IRLS()
    # For PQL fits newsXaug has not been needed in the damping loop but will be needed after exiting this fn
    #   (e.g., for its next call -> LevMarblob <- get_from_MME(sXaug=sXaug, which="LevMar_step", LMrhs=zInfo$scaled_grad, damping=damping))
    if (GLGLLM_const_w) {
      APHLs_args$sXaug <- newsXaug <- sXaug
    } else {
      newsXaug <- do.call(corr_method, mMatrix_arglist)
      if (Trace) { 
        tracechar <- ifelse(.BLOB(newsXaug)$nonSPD,"!",".")
        if (processed$p_v_obj=="p_v") { # v estimation within HL11
          cat(stylefn_v(tracechar))
        } else  cat(stylefn(tracechar)) # PQL/L, vb extimation
      }
      APHLs_args$sXaug <- newsXaug
    } 
    APHLs_args$which <- processed$p_v_obj # "p_v" # 
    newAPHLs <- do.call(".calc_APHLs_from_ZX", APHLs_args) 
  }
  #names(damping) <- which_LevMar_step
  RESU <- list(lik=newlik, APHLs=newAPHLs, damping=damping, 
               sXaug=newsXaug, fitted=fitted, 
               eta=newmuetablob$sane_eta, muetablob=newmuetablob, wranefblob=newwranefblob,
               breakcond=breakcond, 
               v_h=v_h, u_h=u_h, w.resid=neww.resid, weight_X=newweight_X)
  # if ( ! first_it) { # if not break in first iteration
  #   RESU$conv_logL_not_first_it <- conv_logL
  # }
  if ( ! GLMMbool ) {
    RESU$ZAL_scaling <- newZAL_scaling
    RESU$Xscal <- newXscal ## newXscal contains ZAL with new scaling, but (as for any Xscal) independent from weight_X since weight_X is applied only locally in the corr_method
    Vscaled_beta[seq_n_u_h] <- v_h/newZAL_scaling ## represent solution in new scaling...
  } 
  RESU$Vscaled_beta <- Vscaled_beta 
  return(RESU)
}

#copies to allow independent debug()ing
.do_damped_WLS_v_in_b <- .do_damped_WLS 
.do_damped_WLS_outer <- .do_damped_WLS


.WLS_substitute <- function(Xscal, Vscaled_beta, ypos, off, etaFix, seq_n_u_h, ZAL_scaling, mod_attr, 
                            lambda_est, H_global_scale, ZAL, 
                            which_i_llblock, n_u_h, nobs, processed, phi_est, corr_method=processed$corr_method,
                            wranefblob, weight_X, Trace,stylefn) {
  
  # Vscaled_beta must have been provided by something else than damped_WLS_blob
  # drop, not as.vector(): names are then those of (final) eta and mu -> used by predict() when no new data
  fitted <- drop(Xscal %*% Vscaled_beta) ## length nobs+nr ! 
  eta <- fitted[ypos] + off
  RESU <- list()
  if (is.null(etaFix$v_h)) { 
    v_h <- Vscaled_beta[seq_n_u_h] * ZAL_scaling
    if (mod_attr$GLMMbool) {
      RESU$u_h <- RESU$v_h <- v_h ## keep input wranefblob since lambda_est not changed
    } else {
      RESU$u_h <- u_h <- processed$u_h_v_h_from_v_h(v_h)
      if ( ! is.null(maybe <- attr(u_h,"v_h"))) v_h <- maybe
      RESU$v_h <- v_h
      ## update functions u_h,v_h
      RESU$wranefblob <- wranefblob <- processed$updateW_ranefS(u_h=u_h,v_h=v_h, lambda=lambda_est)
      #if ( ! mod_attr$GLMMbool) { # updates ZAL_scaling and functions of it
        RESU$ZAL_scaling <- ZAL_scaling <- 1/sqrt(wranefblob$w.ranef*H_global_scale) ## Q^{-1/2}/s
        RESU$Xscal <- Xscal <- .calc_Xscal_newscaled(Xscal, ZAL_scaling, ZAL, which_i_llblock, 
                                       n_u_h, seq_n_u_h, processed)
        Vscaled_beta[seq_n_u_h] <- v_h/ZAL_scaling ## represent solution in new scaling...
        RESU$Vscaled_beta <- Vscaled_beta 
      #} # Xscal immediately neededfor  sXaug <- do.call(corr_method, list(Xaug=Xscal...
    }
  }
  RESU$muetablob <- muetablob <- .muetafn(eta=eta,BinomialDen=processed$BinomialDen,processed=processed, phi_est=phi_est) 
  fitted[ypos] <- muetablob$sane_eta
  RESU$fitted <- fitted
  if ( (! mod_attr$LLM_const_w) && (! mod_attr$GLGLLM_const_w) ) {
    ## weight_X and Xscal vary within loop if ! LMM since at least the GLMweights in w.resid change
    RESU$w.resid <- .calc_w_resid(muetablob$GLMweights,phi_est, obsInfo=processed$how$obsInfo)
    # at this point w.resid is always the result of .calc_w_resid()
    # and when it is a list with info about mv model it has a complete vector $w_resid.
    H_w.resid <- .calc_H_w.resid(RESU$w.resid, muetablob=muetablob, processed=processed) # for LLF w.resid is not generally defined.
    RESU$weight_X <- .calc_weight_X(Hobs_w.resid=H_w.resid, H_global_scale=H_global_scale, obsInfo=processed$how$obsInfo) ## sqrt(s^2 W.resid)  
    RESU$sXaug <- do.call(corr_method, 
                     list(Xaug=Xscal, weight_X=RESU$weight_X, w.ranef=wranefblob$w.ranef, H_global_scale=H_global_scale,
                          force_QRP=TRUE))
    if (Trace)  if(.BLOB(RESU$sXaug)$nonSPD) {cat(stylefn("!"))} else cat(stylefn("."))
  } ## ergo sXaug is not updated for LMM (no need to)
  return(RESU) ## contains only updated quantities
}

.calc_m_grad_obj <- function(zInfo, 
                             z1_eta=zInfo$z1-etamo, 
                             z1_sscaled_eta=zInfo$z1_sscaled-etamo, # This distinct arg is needed bc I cannot assume the $z1_sscaled is in zInfo (when .calc_m_grad_obj called in spprec LLM)
                             GLMMbool, v_h, wranefblob, 
                             H_w.resid, # z1 (->z1_eta and z1_sscaled_eta) should also have been defined from such weights 
                             # (these are not the gradient weights, BUT cf dlogcLdeta/H_w.resid)
                             #  
                             ZAL, X.pv, etamo) {
  #z1_eta <- zInfo$z1-etamo
  #z1_sscaled_eta <- zInfo$z1_sscaled - etamo # zAug[-seq_n_u_h]-etamo # z_1-sscaled-etamo
  if (GLMMbool) {
    dlogfvdv <-  - v_h * wranefblob$w.ranef
  } else dlogfvdv <- (zInfo$z2 - v_h) * wranefblob$w.ranef
  ## the gradient for -p_v (or -h), independent of the scaling
  # H_w.resid here: cf 'Evaluation of the gradient: total expression'
  m_grad_obj <- c( ## drop() avoids c(Matrix..) 
    m_grad_v <- drop(.crossprod(ZAL, H_w.resid * z1_eta) + dlogfvdv), # Z'W(z_1-eta)+ dlogfvdv 
    drop(.crossprod(X.pv, H_w.resid * z1_sscaled_eta)) # X'W(z_1-sscaled-eta)
  )
  # This can be written
  # drop(crossprod(diag(c(rep(1,n_u_h),H_w.resid)) %*% Xaug,c(drop(.crossprod(ZAL, zInfo$sscaled*H_w.resid)) + dlogfvdv,z1_sscaled_eta)))
  m_grad_obj
}

# code deduced in doc_code/rhs_m_...
.calc_rhs_m_grad_as_szAug <- function(zInfo, 
                                      z1_sscaled_eta, # This distinct arg is needed bc I cannot assume the $z1_sscaled is in zInfo (when .calc_m_grad_obj called in spprec LLM)
                                      z1_eta,
                                      GLMMbool, v_h, wranefblob, 
                                      H_w.resid, 
                                      weight_X, # this comes from sXaug being (diag(c(1/ZAL_scaling,attr(sXaug,"weight_X"))) %*% Xaug %*% diag(c(ZAL_scaling,1)))
                                      H_global_scale,
                                      ZAL_scaling,
                                      ZAL)   {
  
  if (GLMMbool) {
    dlogfvdv <-  - v_h * wranefblob$w.ranef
  } else dlogfvdv <- (zInfo$z2 - v_h) * wranefblob$w.ranef
  rhs_p_v <- sign(H_w.resid)*sqrt(abs(H_w.resid)/H_global_scale) *z1_sscaled_eta
  c( 
    (drop(.crossprod(ZAL, z1_eta*H_w.resid)) + dlogfvdv)  -drop(.crossprod(ZAL, weight_X*rhs_p_v))* ZAL_scaling,
    rhs_p_v
  )
}

# cf doc_code/clean_explanation_scaled_rhs_and_QRP_solution.R for what the resulting scaled_rhs should achieve
# fato_of_gradient.nb was used to get the right algebra here.
.calc_rhs_m_SCALED_grad_as_szAug <- function(zInfo, 
                                             etamo,
                                             z1_sscaled_eta=zInfo$z1_sscaled-etamo, # This distinct arg is needed bc I cannot assume the $z1_sscaled is in zInfo (when .calc_m_grad_obj called in spprec LLM)
                                             z1_eta=zInfo$z1-etamo,
                                             GLMMbool, v_h, wranefblob, 
                                             H_w.resid, 
                                             weight_X, # this comes from sXaug being (diag(c(1/ZAL_scaling,attr(sXaug,"weight_X"))) %*% Xaug %*% diag(c(ZAL_scaling,1)))
                                             H_global_scale,
                                             ZAL_scaling,
                                             ZAL)   {
  
  if (GLMMbool) {
    dlogfvdv <-  - v_h * wranefblob$w.ranef
  } else dlogfvdv <- (zInfo$z2 - v_h) * wranefblob$w.ranef
  if (FALSE) {
    rhs_p_v <- sqrt(abs(H_w.resid)/H_global_scale) *z1_sscaled_eta
    rhs_h_in_m_grad_as_szAug <- (drop(.crossprod(ZAL, z1_eta*H_w.resid)) + dlogfvdv)  -drop(.crossprod(ZAL, weight_X*rhs_p_v))* ZAL_scaling
    H_global_scale*c(
      ZAL_scaling * (rhs_h_in_m_grad_as_szAug + (ZAL_scaling-1)* drop(.crossprod(ZAL, weight_X*rhs_p_v))),
      rhs_p_v
    )
  } else {
    # reduces to
    H_global_scale*c( 
      ZAL_scaling * (drop(.crossprod(ZAL, zInfo$sscaled*H_w.resid)) + dlogfvdv),
      sqrt(abs(H_w.resid)/H_global_scale) *z1_sscaled_eta
    )
  }
}



.solve_IRLS_as_ZX <- 
  function(X.pv=processed$AUGI0_ZX$X.pv, 
           ZAL, y, ## could be taken from processed ? 
           n_u_h=length(u_h), 
           H_global_scale, # from the initial w.resid: not updated when w.resid is locally updated
           lambda_est, muetablob, off=processed$off, maxit.mean, etaFix,
           wranefblob, processed,
           ## supplement for ! LMM
           phi_est, 
           ## supplement for LevenbergM or ! GLMM   ........ and GLMMs ! .makeCovEst1 -> objfn -> .solve_IRLS_as_ZX !
           w.resid,  # distinctly used for the gradient
           H_w.resid,  # distinctly used for the Hessian
           ## supplement for LevenbergM
           beta_eta,
           ## supplement for ! GLMM
           u_h, v_h, # for_init_z_args, 
           ## supplement for intervals
           for_intervals,
           verbose=processed$verbose,
           LevM_HL11_method=.spaMM.data$options$LevM_HL11_method,
           ## ignored but for consistency of arguments with .solve_IRLS_as_spprec:
           corrPars
  ) {
    trace <- verbose["TRACE"]
    if (trace) {
      cat(">") 
      if (verbose["trace"]) cat(.pretty_summ_lambda(lambda_est,processed))
    }
  pforpv <- ncol(X.pv)
  nobs <- length(y)
  seq_n_u_h <- seq_len(n_u_h)
  ypos <- n_u_h+seq_len(nobs)
  lcrandfamfam <- attr(processed$rand.families,"lcrandfamfam")
  mod_attr <- attributes(processed[["models"]])
  LMMbool <-  mod_attr$LMMbool
  GLMMbool <- mod_attr$GLMMbool
  LevenbergM <- (processed$LevenbergM["LM_start"] && is.null(for_intervals))
  is_HL1_1 <- (processed$HL[1L]==1L)
  fpot_cond <- processed$spaMM_tol$fpot_cond
  if (fpot_cond) fpot_tol <- processed$spaMM_tol$fpot_tol
  if ( is.null(for_intervals) && is_HL1_1) {
    if (pforpv==0L) { # outer beta
      which_LevMar_step <- "v" # appeared with outer beta estim
      v_iter <- 0L
    } else which_LevMar_step <- default_b_step <- LevM_HL11_method[["b_step"]] 
    rescue_thr <- processed$spaMM_tol$rescue_thr
    rescue_nbr <- 0L
    prev1_rescued <- FALSE
  } else which_LevMar_step <- "v_b" 
  old_relV_beta <- NULL
  not_moving <- FALSE
  damped_WLS_blob <- NULL
  d_relV_b_tol <- processed$spaMM_tol$d_relV_b_tol 
  d_relV_b_tol_LM <- processed$spaMM_tol$d_relV_b_tol_LM 
  if ( LevenbergM) { 
    dampings_env <- list2env(.spaMM.data$options$spaMM_tol$dampings_env_v)
  } 
  if ( ! LMMbool) {
    checkpot_min_it <- as.integer(maxit.mean/4L) # (see ref to pot4improv in test-mv-nested for a test)
    constant_zAug_args <- list(n_u_h=n_u_h, nobs=nobs, pforpv=pforpv, y=y, off=off, ZAL=ZAL, processed=processed)
    # if ( ! GLMMbool) {
    #   constant_init_z_args <- c(list(lcrandfamfam=lcrandfamfam, nobs=nobs, lambda_est=lambda_est, ZAL=ZAL),  
    #                             # fit_as_ZX args specific for ! GLMM:
    #                             for_init_z_args,
    #                             #
    #                             mget(c("cum_n_u_h","rand.families"),envir=processed))
    # } 
  } 
  
  ##### initial sXaug
  H_global_scale <- .calc_H_global_scale(H_w.resid)
  ZAL_scaling <- 1/sqrt(wranefblob$w.ranef*H_global_scale) ## Q^{-1/2}/s
  Xscal <- .make_Xscal(ZAL, ZAL_scaling = ZAL_scaling, processed=processed, as_matrix=.eval_as_mat_arg(processed))
  which_i_llblock <- .which_i_llblock(Xscal, n_u_h) # preprocessing for faster updating of (sparse) Xscal when scaling changes
  weight_X <- .calc_weight_X(Hobs_w.resid=H_w.resid, H_global_scale=H_global_scale, obsInfo=processed$how$obsInfo) ## sqrt(s^2 [H_]W.resid) # -> .... sqrt([H_]w.resid * H_global_scale)
  corr_method_fn <- get(processed$corr_method,asNamespace("spaMM"), inherits=FALSE)
  sXaug <- corr_method_fn(Xaug=Xscal, weight_X=weight_X, w.ranef=wranefblob$w.ranef, H_global_scale=H_global_scale,
                             force_QRP= ! LevenbergM)
  if (trace) {
    stylefn <- switch(which_LevMar_step,
                      v=.spaMM.data$options$stylefns$vloop,
                      V_IN_B=.spaMM.data$options$stylefns$v_in_loop,
                      .spaMM.data$options$stylefns$betaloop )
    if (LevenbergM) cat("LM")
    if (.BLOB(sXaug)$nonSPD) {cat(stylefn("!"))} else cat(stylefn("."))
  }
  if ( ! is.null(for_intervals)) {
    Vscaled_beta <- c(v_h/ZAL_scaling ,for_intervals$beta_eta)
    fixefobjfn <- names(for_intervals$fixeflik)
  } else if (LevenbergM) {
    Vscaled_beta <- c(v_h/ZAL_scaling ,beta_eta)
  } 
  # to be evaluated once when it becomes needed:
  delayedAssign("constant_v_infer_args", list( # ultimately for the .solve_v_h_IRLS() call
    X.pv=X.pv, ZAL=ZAL, y=y, n_u_h=n_u_h, # H_global_scale=H_global_scale,
    lambda_est=lambda_est, off=off,
    maxit.mean=maxit.mean, # i.e. maxit.mean affects also .solve_v_h_IRLS() calls
    etaFix=etaFix,
    processed=processed, phi_est=phi_est, # for_init_z_args=for_init_z_args,
    trace=trace, dampings_env=dampings_env))
  ## Loop controls:
  allow_LM_restart <- ( ! LMMbool && ! LevenbergM && is.null(for_intervals) && is.na(processed$LevenbergM["user_LM"]) )
  LMcond <- - 10. 
  if (allow_LM_restart) {
    keep_init <- new.env() # ___F I X M E___ that remains an ugly bit of code...
    names_keep <- c("sXaug","wranefblob","muetablob","u_h","w.resid","H_w.resid","v_h","ZAL_scaling","weight_X","Xscal","beta_eta",
                    "old_relV_beta")
    for (st in names_keep) keep_init[[st]] <- environment()[[st]]
  }
  best_HL1_lik <- -Inf
  pot4improv <- NULL
  ################ L O O P ##############
  for (innerj in 1:maxit.mean) {
    if( ! LevenbergM && allow_LM_restart) { ## FIXME the next step improvement would be 
      #  ./. to keep track of lowest lambda that created problem and use LM by default then
      # if (innerj>1L && notSPD) {
      #   LevenbergM <- TRUE
      #   for (st in names_keep) assign(st,keep_init[[st]])
      #   Vscaled_beta <- c(v_h/ZAL_scaling ,beta_eta) ## RHS from assign(st,keep_init[[st]])
      #   dampings_env <- list2env(.spaMM.data$options$spaMM_tol$dampings_env_v)
      #   damped_WLS_blob <- NULL
      #   allow_LM_restart <- FALSE 
      #   if ( which_LevMar_step=="v_b") { 
      #     ## The LevM.negbin test finds "strict_v|b" poorer than "V_IN_B"" (note some divergent p_v's)  which led to:
      #     # which_LevMar_step <- "V_IN_B" ## not modified by if (... ! is.null(damped_WLS_blob) ...) before being used.
      #     # BUT the  optim_LevM's (update(br$fullfit,fixed=... test shows one should keep using "v_b" here.
      #     # otherwise "!LM" differs (and is poorer) from LevM=TRUE (which indeed starts from "v_b")
      #     # However, it's not clear why "V_IN_B" is poorer (and it's not the step on which the loop terminates)  => _F I X M E_ 
      #     strictv_parent_info <- c(from=which_LevMar_step, breakcond="")
      #     rescue_thr <- processed$spaMM_tol$rescue_thr
      #     rescue_nbr <- 0L
      #     prev1_rescued <- FALSE
      #   }
      # } else 
      if (innerj>3L) {
        crit <- abs_d_relV_beta/(old_abs_d_relV_beta+1e-8)
        LMcond <- LMcond + mean(sqrt(crit))^2 
        #previously, mean(crit) could be upset by *one* small old_abs_d_relV_beta() value ()eg glmmTMB COMPoisson example)
        # conversely, median is too insensitive (long tests globally)
        # we could adapt the "norm()" to the expected costs/benefit of different thresholds (higher cost for COMPoisson) __F I X M E__ 
        if (LMcond/innerj>0.5 || 
            (innerj> checkpot_min_it && pot4improv > max(10,old_pot4improv))
            ) {
          LevenbergM <- TRUE
          if (trace) cat("!LM") 
          for (st in names_keep) assign(st,keep_init[[st]])
          Vscaled_beta <- c(v_h/ZAL_scaling ,beta_eta) ## RHS from assign(st,keep_init[[st]])
          dampings_env <- list2env(.spaMM.data$options$spaMM_tol$dampings_env_v)
          damped_WLS_blob <- NULL
          allow_LM_restart <- FALSE 
          if ( which_LevMar_step=="v_b") { 
            ## The LevM.negbin test finds "strict_v|b" poorer than "V_IN_B"" (note some divergent p_v's)  which led to:
            # which_LevMar_step <- "V_IN_B" ## not modified by if (... ! is.null(damped_WLS_blob) ...) before being used.
            # BUT the  optim_LevM's (update(br$fullfit,fixed=... test shows one should keep using "v_b" here.
            # otherwise "!LM" differs (and is poorer) from LevM=TRUE (which indeed starts from "v_b")
            # However, it's not clear why "V_IN_B" is poorer (and it's not the step on which the loop terminates)  => _F I X M E_ 
            strictv_parent_info <- c(from=which_LevMar_step, breakcond="")
            rescue_thr <- processed$spaMM_tol$rescue_thr
            rescue_nbr <- 0L
            prev1_rescued <- FALSE
          }
        }
      }
      if (innerj>2L) {
        old_abs_d_relV_beta <- abs_d_relV_beta
        old_pot4improv <- pot4improv
      }
    }
    ##### get the lik of the current state
    if ( ! is.null(for_intervals)) {
      loc_logLik_args <- list(sXaug=sXaug, processed=processed, phi_est=phi_est,
                              lambda_est=lambda_est, dvdu=wranefblob$dvdu, u_h=u_h, muetablob=muetablob)
      oldlik <- unlist(do.call(".calc_APHLs_from_ZX",loc_logLik_args)[fixefobjfn]) # unlist keeps name
    } else if (LevenbergM) { ## then logL is necessary to check for increase
      if (is.null(damped_WLS_blob)) {
        oldAPHLs <- .calc_APHLs_from_ZX(sXaug=sXaug, processed=processed, phi_est=phi_est, which=processed$p_v_obj, 
                                        lambda_est=lambda_est, dvdu=wranefblob$dvdu, u_h=u_h, 
                                        muetablob=muetablob)
      } else { ## Levenberg and innerj>1
        oldAPHLs <- damped_WLS_blob$APHLs
      }
    } 
    #####
    
    ##### RHS
    if (LMMbool) {
      wzAug <- c(rep(0,n_u_h),(y-off)*weight_X)
    } else {
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
                               sXaug=sXaug, 
                               w.ranef=wranefblob$w.ranef, 
                               w.resid=w.resid,
                               ############################# ZAL_scaling=ZAL_scaling,
                               z2=z2) )
      zInfo <- do.call(".calc_zAug_not_LMM",calc_zAug_args) ## dlogfvdv is represented indirectly in $z2, as shown by code in .calc_m_grad_obj()
      wzAug <- c(zInfo$y2_sscaled/ZAL_scaling, (zInfo$z1_sscaled)*weight_X) 
    }
    ## keep name 'w'zAug to emphasize the distinct weightings  of zaug and Xaug (should have been so everywhere)
    #####
    ##### improved  Vscaled_beta   
    if ( ! is.null(for_intervals)) {
      currentDy <- (for_intervals$fixeflik-oldlik)
      # given current focal-parameter estimate and associated current lik, 
      # try to guess the new estimate that will bring the likelihood closer to the CI threshold:
      intervalBlob <- .intervalStep_ZX(old_Vscaled_beta=Vscaled_beta,
                                       sXaug=sXaug,szAug=wzAug,
                                       for_intervals=for_intervals,
                                       currentlik=oldlik,currentDy=currentDy) # return value is list(Vscaled_beta=...)
      damped_WLS_blob <- NULL
      Vscaled_beta <- intervalBlob$Vscaled_beta
    } else if (LevenbergM) { ## excludes IRLS
      ## (w)zAug is all what is needed for the direct solution of the extended system. in GLMM case
      # Hence wZaug contains Phi z_2 including (Phi v^0 +dlogfvdv)/ZAL_scaling (from components of hlik)
      ## now we want the LHS of a d_beta_v solution
      m_grad_obj <- .calc_m_grad_obj(zInfo, GLMMbool=GLMMbool, v_h=v_h, wranefblob=wranefblob, 
                                     H_w.resid=.BLOB(sXaug)$H_w.resid, 
                                     ZAL=ZAL, X.pv=X.pv,etamo=muetablob$sane_eta - off)
      if (trace>1L) {
        stylefn <- switch(which_LevMar_step,
                       v=.spaMM.data$options$stylefns$vloop,
                       V_IN_B=.spaMM.data$options$stylefns$v_in_loop,
                       .spaMM.data$options$stylefns$betaloop )
        if (pforpv) { 
          maxs_grad <- c(max(abs(m_grad_obj[seq_n_u_h])),max(abs(m_grad_obj[-seq_n_u_h])))
        } else maxs_grad <- c(max(abs(m_grad_obj[seq_n_u_h])), 0) #outer beta
        # if (innerj==44L) browser()
        cat(stylefn("iter=",innerj,", max(|grad|): v=",maxs_grad[1L],"beta=",maxs_grad[2L],";\n"))
      }
      zInfo$gainratio_grad <- m_grad_obj ## before rescaling
      # gradient for scaled system from gradient of objective
      scaled_grad <- H_global_scale * m_grad_obj
      scaled_grad[seq_n_u_h] <- scaled_grad[seq_n_u_h] * ZAL_scaling 
      zInfo$scaled_grad <- scaled_grad
      constant_APHLs_args <- list(processed=processed, which=processed$p_v_obj, sXaug=sXaug, phi_est=phi_est, lambda_est=lambda_est)
      # the following block needs m_grad_v the new m_grad_v hence its position
      if (is_HL1_1 && ! is.null(damped_WLS_blob)) {
        #### Get next LM step && conditionally update old_relV_beta ####
        # we assess convergence at the end of the loop by comparing old_relV_beta to relV_beta. We update old_relV_beta 
        #    (1) here, in all cases where v has been updated;
        # or (2) in one reversal case, it needs to be updated after the test at the end of the loop.
        if (just_rescued <- identical(attr(damped_WLS_blob, "step"), "rescue")) {
          rescue_nbr <- rescue_nbr + 1L
          old_relV_beta <- relV_beta 
          if (pforpv) { 
            if (prev1_rescued || rescue_nbr > rescue_thr["V_IN_B"]) {
              which_LevMar_step <- "V_IN_B" # # [yellow [cyan .cyan .....ul green...
            } else which_LevMar_step <- "v_b" # outer beta
          }
        } else if (which_LevMar_step=="v_b" || which_LevMar_step=="b_from_v_b" ) { 
          if (rescue_nbr > rescue_thr["strictv"]  &&  #rescue has been previously needed in the outer loop
              damped_WLS_blob$breakcond != "low_pot" ## in particular, if OK_gain, we play safer if we know a problem occurred previously
                                                      ## "stuck_obj" occurs here too...
          ) { 
            strictv_parent_info <- c(from=which_LevMar_step, breakcond=damped_WLS_blob$breakcond)
            which_LevMar_step <- "strict_v|b" # 
            # } else if (max(abs(m_grad_v)) > max(abs(old_m_grad_v))) which_LevMar_step <- "v" # test is fausse bonne idee...
          } else {
            v_parent_info <- c(from=which_LevMar_step, breakcond=damped_WLS_blob$breakcond)
            # Tried this from [v3.12.34 -> 54] with poor effect on test-nloptr #362... remains poor when tested from v3.13.11
            # if (damped_WLS_blob$breakcond=="stuck_obj") {
            #   which_LevMar_step <- "V_IN_B"
            # } else
               which_LevMar_step <- "v" # standard switch from yellow to underlined cyan. It is generally not a good idea to switch immediately to "V_IN_B" 
            v_iter <- 0L
          }
        } else if (which_LevMar_step=="v") {
          if (damped_WLS_blob$breakcond == "low_pot") { ## LevMar apparently maximized h wrt v after several iterations
            #cat(damped_WLS_blob$breakcond)
            old_relV_beta <- relV_beta ## serves to assess convergence !!! which is thus dependent on condition ( hlik_stuck || ! need_v_step)
            if (pforpv) which_LevMar_step <- default_b_step # We should not reach this line when RHS is "v_in_b"  # ! outer beta
          } else {
            v_iter <- v_iter+1L
            # Series of "v" are series of small gains ("OK_gain"). seen in test_negbin_difficult, not in test nloptr nor COMPoisson.
            # Assuming a geometric gain we could estimate their time until lowpot is reached
            #  and determine from it the condition for switching: (tried: this fails massively)
            if (v_iter>10L && pforpv) which_LevMar_step <- "V_IN_B" # [yellow [cyan .cyan .....ul green... Ie a [[v_h_IRLS]+probable beta updating]
          }
        } else if (which_LevMar_step=="strict_v|b") {
          old_relV_beta <- relV_beta 
          if (strictv_parent_info[["from"]] %in% c("v_b","b_from_v_b" )
              && strictv_parent_info[["breakcond"]] !="low_pot") { # (which implies that rescue was not just called)
            which_LevMar_step <- "v_b"
          } else { # strictv was called after V_IN_B (!)
            if (rescue_nbr > rescue_thr["re_V_IN_B"]) {
              which_LevMar_step <- "V_IN_B"
            } else which_LevMar_step <- "v_b"  
          }
        } else if (which_LevMar_step=="V_IN_B") { 
          ## .wrap_do_damped_WLS_outer(which_LevMar_step=="V_IN_B") has run .do_damped_WLS_outer(which_LevMar_step=="b_&_v_in_b") ie.
          ## a damping loop of {a beta updating followed by a v_h_IRLS}'s. We assess the resulting breakcond.
          breakcond <- damped_WLS_blob$breakcond
          if (breakcond=="stuck_obj" || breakcond=="div_gain") {
            strictv_parent_info <- c(from=which_LevMar_step,breakcond=breakcond)
            which_LevMar_step <- "strict_v|b" # the call to "strict_v|b" may seem odd but results in clean optim
            #If we did that in the wrap... then we would next compare two identical "strict_v|b" 
          } else { # a good breakcond MAY result in switching back to v_b
            old_relV_beta <- relV_beta 
            if (rescue_nbr > rescue_thr["re_V_IN_B"]) {
              which_LevMar_step <- "V_IN_B"
            } else which_LevMar_step <- "v_b" 
          }
        } else if (default_b_step=="v_in_b") { # presumably not used
          old_relV_beta <- relV_beta
        } else { ## "b" or any unanticipated case # presumably not used 
          warning("Unexpected case in .solve_IRLS_as_ZX(): please contact the package maintainer.")
          # as v_b case:
          if (rescue_nbr > rescue_thr["strictv"] &&  #rescue has been previously needed in the outer loop
              damped_WLS_blob$breakcond != "low_pot") { 
            strictv_parent_info <- c(from=which_LevMar_step,breakcond=damped_WLS_blob$breakcond)
            which_LevMar_step <- "strict_v|b" 
          } else {
            v_parent_info <- c(from=which_LevMar_step, breakcond=damped_WLS_blob$breakcond)
            which_LevMar_step <- "v" 
            v_iter <- 0L
          }
        }
        prev1_rescued <- just_rescued 
      } 
      new_damping <- .get_new_damping(dampings_env$v[[which_LevMar_step]], which_LevMar_step)
      damped_WLS_blob <- .wrap_do_damped_WLS_outer(
        damped_WLS_fn = .do_damped_WLS_outer,
        LevM_HL11_method=LevM_HL11_method,
        rescue= (is_HL1_1 && rescue_thr["rescue"]),
        which_LevMar_step=which_LevMar_step,
        old_relV_beta=old_relV_beta,
        sXaug=sXaug, zInfo=zInfo, 
        old_Vscaled_beta=Vscaled_beta,
        oldAPHLs=oldAPHLs,
        APHLs_args = constant_APHLs_args,
        damping=new_damping,
        Trace= trace,
        ypos=ypos,off=off,
        GLMMbool=GLMMbool,etaFix=etaFix,
        lambda_est=lambda_est,
        wranefblob=wranefblob,seq_n_u_h=seq_n_u_h,ZAL_scaling=ZAL_scaling,
        processed=processed, Xscal=Xscal,
        phi_est=phi_est, H_global_scale=H_global_scale, n_u_h=n_u_h, ZAL=ZAL,
        which_i_llblock=which_i_llblock,
        constant_v_infer_args=constant_v_infer_args,
        looseness= if ( is.null(damped_WLS_blob) ||  ## start strict
                        new_damping>1e-7) {## use strict when there are trace of difficulties (in particular, failure to improve) 
          1 } else {processed$spaMM_tol$loose_fac},
        low_pot=NULL ## explicit for clarity, but its the default
      ) 

      #old_m_grad_v <- m_grad_v
      #damping_barring_rescue <- damped_WLS_blob$damping
      #cat(crayon::yellow(paste0(names(damping_barring_rescue),"=",damping_barring_rescue[[1]])))
      #dampings_env$v[[names(damping_barring_rescue)]] <- damping_barring_rescue[[1]]
      dampings_env$v[[attr(damped_WLS_blob,"step")]] <- damped_WLS_blob$damping 
      # : when step is rescue the rescued step does not provide an updated damping; dampings_env$v[["rescue"]] is ignored. 
      ## LevM PQL
      if (! is_HL1_1) {
        if (damped_WLS_blob$lik < oldAPHLs$hlik) { ## if LevM step failed to find a damping that increases the hlik
          # Tis should occur only bc of (1) numerically challenging conditions e.g mu close to bounds; or (2) optimum has been 
          # found and floating point inaccuracies matter. We try to exclude the second case by the following test: 
          if ( ! ((breakcond <- damped_WLS_blob$breakcond)=="low_pot" && attr(breakcond,"very_low_pot"))) {
            damped_WLS_blob <- NULL
            .diagnose_conv_problem_LevM( beta_cov_info=get_from_MME(sXaug,which="beta_cov_info_from_sXaug"), processed) # writes into 'processed'
            wzAug <- c(zInfo$y2_sscaled/ZAL_scaling, (zInfo$z1_sscaled)*weight_X)
            Vscaled_beta <- get_from_MME(sXaug,szAug=wzAug) # vscaled= v scaling so that v has 'scale' * ZAL_scaling
            if (TRUE) { # not clear what is best here
              break
            } else {
              LevenbergM <- FALSE ## desperate move... 
              # for (st in names_keep) assign(st,keep_init[[st]])
              # Vscaled_beta <- c(v_h/ZAL_scaling ,beta_eta) ## RHS from assign(st,keep_init[[st]])
            }
          }
        } 
      }
    # } else if (.BLOB(sXaug)$nonSPD) { ## IRLS: always accept new v_h_beta
    #   Vscaled_beta <- c(v_h/ZAL_scaling ,beta_eta)
    #   if (trace) cat("!")
    #   damped_WLS_blob <- NULL
    } else {
      ################ direct solution for Vscaled_beta
      if (FALSE && ! is.null(.BLOB(sXaug)$signs)) { # This code is correct (spcorr QRP), 
        # explained in doc_code/clean_explanation_scaled_rhs_and_QRP_solution.R
        # hence get_from_MME(sXaug,szAug=wzAug) should give the same result.
        scaled_rhs <- .calc_rhs_m_SCALED_grad_as_szAug(
          zInfo = zInfo,
          etamo=  muetablob$sane_eta - off,
          #z1_sscaled_eta = zInfo$z1_sscaled-etamo, z1_eta = zInfo$z1-etamo, 
          GLMMbool = GLMMbool,
          v_h = v_h,wranefblob = wranefblob,
          H_w.resid = .BLOB(sXaug)$H_w.resid,
          H_global_scale=H_global_scale, ZAL_scaling=ZAL_scaling,
          weight_X=weight_X, ZAL = ZAL)
        augsigns <- c(rep(1,attr(sXaug,"n_u_h")), .BLOB(sXaug)$signs) 
        if (innerj>1L) {
          Vscaled_beta <- drop(get_from_MME(sXaug,szAug = scaled_rhs + sXaug %*% Vscaled_beta))
        } else Vscaled_beta <- get_from_MME(sXaug,szAug = scaled_rhs + c(1/ZAL_scaling,attr(sXaug,"weight_X"))*c(v_h, # input 'v_h' are the unscaled, userlevel ones.
                                                                                                                 muetablob$sane_eta))
        # V'"scaled_beta"'[seq_n_u_h] <- Vscaled_beta[seq_n_u_h]*ZAL_scaling # To get same values as in spprec.
        # BUT spprec has trivial ZAL_scaling=1 so it's normal that its 'scaled' v_h as distinct. The meaning is actually
        # unscaled_v_h <- Vscaled_beta[seq_n_u_h]*ZAL_scaling # To get same values as in spprec.
      } else Vscaled_beta <- get_from_MME(sXaug,szAug=wzAug) # vscaled= v scaling so that v has 'scale' * ZAL_scaling
      ## solution should be equivalent to
      # drop(solve(crossprod(sXaug), 
      #            ( H_global_scale*m_grad_obj * c(ZAL_scaling,1)+ # <= provides d_Vscaled_beta. Note the scalingS on m_grad_obj
      #                crossprod(sXaug , # <= second term to recover the resident Vscaled_beta (at least if we used eta instead of sane_eta...)
      #                          c(1/ZAL_scaling,weight_X) *drop(c(v_h,muetablob$sane_eta-off)))) ) )
      ## Since the second term is not used in the spprec methods, the methods diverge in case of sanitization.
      ## And although spprec may seem more exact, sanitization may well improve convergence...
      damped_WLS_blob <- NULL
      if ( ! LMMbool && innerj>= checkpot_min_it) {
        m_grad_obj <- .calc_m_grad_obj(zInfo, GLMMbool=GLMMbool, v_h=v_h, wranefblob=wranefblob, 
                                       H_w.resid=.BLOB(sXaug)$H_w.resid, 
                                       ZAL=ZAL, X.pv=X.pv,etamo=muetablob$sane_eta - off)
        pot4improv <- get_from_MME(sXaug=sXaug, which="Mg_solve_g", B=m_grad_obj)
      }
    }
    if (trace>5L) .prompt()
    ##### Everything that is needed for 
    #  (1) assessment of convergence: c(v_h*sqrt(wranefblob$w.ranef),beta_eta)
    #  (2) all return elements are updated as function of the latest Vscaled_beta 
    #                               (itself possible updated to new scaling by the following assign()'s).
    #      In particular We need muetablob and (if ! LMM) sXaug, hence a lot of stuff.
    #  Hence, the following code is useful whether a break occurs or not. 
    if ( is.null(damped_WLS_blob) ) { ## fits nothing, but updates variables in case of standard IRLS, or of intervals
      WLS_blob <- .WLS_substitute(Xscal, Vscaled_beta, ypos, off, etaFix, seq_n_u_h, ZAL_scaling, mod_attr=mod_attr, 
                                  lambda_est=lambda_est, H_global_scale, ZAL, 
                                  which_i_llblock, n_u_h, nobs, processed, phi_est=phi_est, 
                                  wranefblob=wranefblob, weight_X=weight_X, Trace=trace, stylefn=stylefn)
      if (fpot_cond && # FALSE except in possible private usage, 
          ! LMMbool  && is.null(for_intervals)) {  
        old_m_grad_obj <- .calc_m_grad_obj(zInfo, GLMMbool=GLMMbool, 
                                           v_h=v_h, wranefblob=wranefblob, 
                                           H_w.resid=.BLOB(sXaug)$H_w.resid, 
                                           ZAL=ZAL, X.pv=X.pv, etamo =muetablob$sane_eta - off)
        Mg_solve_g <- get_from_MME(sXaug=sXaug, which="Mg_solve_g", B=old_m_grad_obj) 
        for (st in names(WLS_blob)) assign(st,WLS_blob[[st]]) #                                     ASSIGN. Includes (at least for GLMMs) w.resid AND H_w.resid
        if (250*Mg_solve_g < fpot_tol) break
      } else  list2env(WLS_blob, envir = environment()) # for (st in names(WLS_blob)) assign(st,WLS_blob[[st]]) 
    } else {
      if (is_HL1_1 && (
        (which_LevMar_step =="V_IN_B" && damped_WLS_blob$breakcond=="OK_gain") ||
        (which_LevMar_step =="v" && damped_WLS_blob$breakcond=="low_pot")
      ) && damped_WLS_blob$APHLs$p_v>best_HL1_lik ) {
        #cat(crayon::red("ICI!\n"))
        best_HL1_damped_WLS_blob <- damped_WLS_blob
        best_HL1_lik <- best_HL1_damped_WLS_blob$APHLs$p_v
        #print(best_HL1_lik)
      } 
      varnames <- intersect(names(damped_WLS_blob), # sXaug (and the weights) need not be present if(GLGLLM_const_w) 
                            c("w.resid", ## !important! cf test-adjacency-corrMatrix.R
                              "weight_X", "Vscaled_beta","wranefblob","v_h","u_h","muetablob", "sXaug"))
      list2env(damped_WLS_blob[varnames], envir = environment()) 
      # for (st in intersect(names(damped_WLS_blob), # sXaug (and the weights) need not be present if(GLGLLM_const_w) 
      #                      c("w.resid", ## !important! cf test-adjacency-corrMatrix.R
      #                        "weight_X", 
      #                        "Vscaled_beta","wranefblob","v_h","u_h","muetablob",
      #                        "sXaug"))) assign(st,damped_WLS_blob[[st]])
      if ( ! GLMMbool ) {
        Xscal <- damped_WLS_blob$Xscal ## contains ZAL with new scaling, but weight_X is not applied since it is applied only locally in the corr_method
        ZAL_scaling <- damped_WLS_blob$ZAL_scaling
      }
    }
    #  At this point all return elements are updated as function of the latest Vscaled_beta.
    #  In particular We need muetablob and (if ! LMM) sXaug, hence a lot of stuff.
    #####
    beta_eta <- Vscaled_beta[n_u_h+seq_len(pforpv)]
    ##### assessment of convergence
    if (innerj==maxit.mean) {
      if (maxit.mean>1L) { 
        if (LevenbergM) {
          processed$LevenbergM["LM_start"] <- TRUE # _F I X M E_ think again
        } ## but by default the local 'LevenbergM' is FALSE and will not become true until divergence is detected
        # print(c(which_LevMar_step,damped_WLS_blob$breakcond))
        if (trace) {
          cat(crayon::red("!"))
        } else if ( ! identical(processed$warned_maxit_mean, TRUE)) {
          processed$warned_maxit_mean <- TRUE
          if (!is.null(for_intervals)) {
            message("Iterative algorithm converges slowly.")
          } else message("Iterative algorithm converges slowly. See help('convergence') for suggestions.")
        }
      }
      break
    } else {    
      # if (innerj==maxit.mean-1L)  cat(crayon::yellow(c(which_LevMar_step,damped_WLS_blob$breakcond))) # browser()
      relV_beta <- c(v_h*sqrt(wranefblob$w.ranef),beta_eta)  ## convergence on v_h relative to sqrt(w.ranef)
      abs_d_relV_beta <- abs(relV_beta - old_relV_beta) ## for ML, comparison between estimates when ( hlik_stuck || ! need_v_step )
      ## abs_d_relV_beta is needed outside this block, and old_relV_beta updated after v_h updating.
      ## Attempt: strongly degrades the perf of nloptr/bobyqa in nloptr test 
      # dlogLdvb <- crossprod(sXaug,wzAug) # is that correct ?
      # abs_d_relV_beta <- drop(100*abs_d_relV_beta*dlogLdvb) # that would actually be d_logL... OK we may need to sum() it... (__F I X M E__)
      if (is_HL1_1 && LevenbergM) {
        # if((which_LevMar_step %in% c("v", "V_IN_B","strict_v|b") || default_b_step=="v_in_b") &&
        #    ! is.null(old_relV_beta)) {
        #   cat(abs(damped_WLS_blob$APHLs[[processed$p_v_obj]])," ")
        # }
        not_moving <- (
          # exclude cases of possible p_v overfit by v_h, such as cases "b" "v_b" "b_from_v_b" "b_&_v_in_b"
          (which_LevMar_step %in% c("v", "V_IN_B","strict_v|b") || default_b_step=="v_in_b") &&
            ( ! is.null(old_relV_beta)) && 
            {
              #print(c(mean(abs_d_relV_beta[seq_n_u_h]),relV_beta[-seq_n_u_h],old_relV_beta[-seq_n_u_h]))
              meanmean <- mean(c(mean(abs_d_relV_beta[seq_n_u_h]), mean(abs_d_relV_beta[-seq_n_u_h])), na.rm=TRUE) # second mean may be NaN
              meanmean < d_relV_b_tol_LM
            }
        )
      } else not_moving <- ( 
        ( ! is.null(old_relV_beta)) && 
        {
          meanmean <- mean(c(mean(abs_d_relV_beta[seq_n_u_h]), mean(abs_d_relV_beta[-seq_n_u_h])), na.rm=TRUE) # second mean may be NaN
          meanmean < d_relV_b_tol
        }
      ) #In ! LevM, v&b are fitted simultaneously without damping
      if (not_moving) {
        # not_moving_Wattr <- .diagnose_coeff_not_moving(coeff_not_moving = not_moving,relV_beta, damped_WLS_blob, innerj, 
        #                                                damping, is_HL1_1, oldAPHLs, Ftol=processed$spaMM_tol$Ftol_LM, trace, LevenbergM,stylefn=identity)
        if (FALSE) {
          etamo <- muetablob$sane_eta - off
          m_grad_obj <- .calc_m_grad_obj(zInfo, z1_eta=zInfo$z1-etamo, z1_sscaled_eta=zInfo$z1_sscaled - etamo, GLMMbool, v_h, wranefblob, 
                                         H_w.resid=.BLOB(sXaug)$H_w.resid, 
                                         ZAL, X.pv)
          if (250*get_from_MME(sXaug=sXaug, which="Mg_solve_g", B=m_grad_obj)<1e-10) break
        } else break
      } # else if (F_I_X_M_E && innerj==maxit.mean-1L) browser()
      # More ad hoc breaks for cases where the coefficients keep moving although the total potential is low:
      if ( ! is.null(damped_WLS_blob) ) {
        if (is_HL1_1) {
          if ( which_LevMar_step=="v_b" && damped_WLS_blob$breakcond=="low_pot" && attr(damped_WLS_blob$breakcond,"no_overfit")) {
            break # motivated by BINARYboot (e.g., replicate 362)
          } else if ( which_LevMar_step=="v" && 
                      damped_WLS_blob$breakcond=="low_pot" &&
                      ( ( ! pforpv ) ||  # outer beta
                        v_parent_info[["breakcond"]]=="low_pot" ## essential condition; otherwise poor fits (eg test-nloptr # 422 p_v_out_def   ) 
                      )
          ) {
            break 
          } else if ( which_LevMar_step=="v_b" && 
                      damped_WLS_blob$breakcond=="stuck_obj" ) { 
            # I removed this break from [v3.12.34 -> 54] with poor effect on test-nloptr #362... remains poor when tested from v3.13.11
            break # this case occurs in the test-nloptr tests
            # not an obvious termination condition, but seems OK following "v"&&"low_pot" (and stops the alternation between these two states)
            # but if previous steps were "V_IN_B" & OK_gain, it may be worth returning to such a step
          } else if ( which_LevMar_step=="strict_v|b" && 
                      # => so that the $breakcond is { "v|b_no_loop" (i.e. no damping loop in this case), with "v_pot4improv" attribute tested next:} 
                      strictv_parent_info[["breakcond"]]=="stuck_obj" && # we don't test "low_pot" here bc this case does not occur according to previous code
                      attr(damped_WLS_blob$breakcond,"v_pot4improv")<1e-10) {
            break # this case occurs in the test-nloptr tests
          } else if ( which_LevMar_step =="V_IN_B" && damped_WLS_blob$breakcond=="low_pot") {
            # I could further test attr(damped_WLS_blob$breakcond,"very_low_pot") here
            break # motivated by poisson 'smaller' fit in test_COMPoisson_difficult with control.HLfit=list(max.iter.mean=1000),
          } #else {print(which_LevMar_step); str(damped_WLS_blob$breakcond)}
        } else {# tendency of PQL/L meanmean to converge very slowly despite low_pot (and indeed tiny changes in hlik) => ad hoc adjustment
          if (damped_WLS_blob$breakcond=="low_pot" && attr(damped_WLS_blob$breakcond,"very_low_pot")) break
        }
      }
      #
      if ( ! (is_HL1_1 && LevenbergM)) { ## This is for the special case of reversal of LevenbergM condition from F to T in  LevM PQL !!!!
        old_relV_beta <- relV_beta
      } ## ELSE old_relV_beta was updated when required in block for which_LevMar_step.
    } 
  } ################ E N D LOOP ##############
  #if (trace>4L) browser() 
  if ( ! is.null(for_intervals) && for_intervals$phi_pred_OK) {
    warnobjfn <- names(for_intervals$warnlik)
    warnlik <- unlist(do.call(".calc_APHLs_from_ZX",loc_logLik_args)[warnobjfn])  
    if ((for_intervals$warnlik-warnlik) < -1e-4 && 
        (is.null(bestlik <- processed$envir$confint_best$lik) || warnlik > bestlik)) {
      if (is.null(bestlik)) {
        locmess <- paste("A higher",warnobjfn,"was found than for the original fit.",
                         "\nThis suggests the original fit did not fully maximize",warnobjfn,
                         "\n (REML fits, or numerical accuracy issues). Expect more information at end of computation.")
        message(locmess)
      }
      processed$envir$confint_best$lik <- warnlik
      processed$envir$confint_best$beta_eta <- .unscale(X.pv, Vscaled_beta[n_u_h+seq_len(pforpv)])
      processed$envir$confint_best$ranPars <- for_intervals$ranFix
      processed$envir$confint_best$ranPars$lambda <- unique(lambda_est) # ugly but get_inits_from_fit() cannot be used in this context...
    }
  }
  if (trace>1L && (LevenbergM))  {
    stylefn <- .spaMM.data$options$stylefns$betalast
    if (pforpv) { # outer beta
      maxs_grad <- c(max(abs(m_grad_obj[seq_n_u_h])),max(abs(m_grad_obj[-seq_n_u_h])))
    } else maxs_grad <- c(max(abs(m_grad_obj[seq_n_u_h])), 0)
    cat(stylefn("iter=",innerj,", max(|grad|): v=",maxs_grad[1L],"beta=",maxs_grad[2L],";\n"))
  }
  if (best_HL1_lik > -Inf) {
    if (best_HL1_lik> damped_WLS_blob$APHLs$p_v) {
      # cat("restoring better fit")
      damped_WLS_blob <- best_HL1_damped_WLS_blob
      varnames <- intersect(names(damped_WLS_blob), # sXaug (and the weights) need not be present if(GLGLLM_const_w) 
                            c("w.resid", ## !important! cf test-adjacency-corrMatrix.R
                              "weight_X", "Vscaled_beta","wranefblob","v_h","u_h","muetablob", "sXaug"))
      list2env(damped_WLS_blob[varnames], envir = environment()) 
    }
  }
  names(beta_eta) <- colnames(X.pv)
  if (! is.null(damped_WLS_blob)) { 
    fitted <- damped_WLS_blob$fitted
    weight_X <- damped_WLS_blob$weight_X ## F I X M E it seems better to store  weight_X  as attr(sXaug,...) and no weight_X elsewhere in output 
  } # otherwise 'fitted' will come from assigning the results of .WLS_substitute()
  RESU <- list(sXaug=sXaug, # updated through assign()
               ## used by calc_APHLs_from_ZX: (in particular can use Vscaled values contained in fitted)
               fitted=fitted, 
               weight_X=weight_X, 
               nobs=nobs, pforpv=pforpv, seq_n_u_h=seq_n_u_h, u_h=u_h, 
               muetablob=muetablob, 
               lambda_est=lambda_est,
               phi_est=phi_est,
               ## used by other code
               beta_eta=beta_eta, w.resid=w.resid, wranefblob=wranefblob, 
               v_h=v_h, eta=muetablob$sane_eta, innerj=innerj
               )
  ## for diagnostic purposes
  if ( ! LMMbool ) {
    RESU$m_grad_obj <- zInfo$gainratio_grad ## frm ZInfo bc other copies of m_grad_obj may be missing if not LevenbergM
  }
  return(RESU)
}

.intervalStep_ZX <- function(old_Vscaled_beta,sXaug,szAug,
                             currentlik, ## not currently used (has been, for diagnostic printing)
                             for_intervals,currentDy) {
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
    targetDy <- (for_intervals$fixeflik-for_intervals$targetlik)
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
  return(list(Vscaled_beta=Vscaled_beta)) # levQ is presumably always dense
}

.intervalStep_glm <- function(old_beta,sXaug,szAug,currentlik,for_intervals,currentDy) {
  ## voir code avant 18/10/2014 pour une implem rustique de VenzonM pour debugage  
  ## somewhat more robust algo (FR->FR: still improvable ?), updates according to a quadratic form of lik near max
  ## then target.dX = (current.dX)*sqrt(target.dY/current.dY) where dX,dY are relative to the ML x and y 
  ## A nice thing of this conception is that if the target lik cannot be approached, 
  ##   the inferred x converges to the ML x => this x won't be recognized as a CI bound (important for locoptim) 
  # if (is.null(names(szAug))) stop("Programming error: 'szAug' must have names")
  if (is.null(colnames(sXaug))) stop("Programming error: 'sXaug' must have colnames")
  parmcol_X <- for_intervals$parmcol_X
  beta <- rep(NA,length(old_beta))
  if (currentDy <0) { 
    beta[parmcol_X] <- old_beta[parmcol_X]
  } else {
    currentDx <- (old_beta[parmcol_X]-for_intervals$MLparm)
    targetDy <- (for_intervals$fixeflik-for_intervals$targetlik)
    Dx <- currentDx*sqrt(targetDy/currentDy)
    ## pb is if Dx=0 , Dx'=0... and Dx=0 can occur while p_v is still far from the target, because other params have not converged.
    ## FR->FR patch:
    if (currentDy<targetDy) { ## we are close to the ML: we extrapolate a bit more confidently
      min_abs_Dx <- for_intervals$asympto_abs_Dparm/1000
    } else min_abs_Dx <- 1e-6 ## far from ML: more cautious move our of Dx=0
    Dx <- sign(currentDx)*max(abs(Dx),min_abs_Dx)
    beta[parmcol_X] <- for_intervals$MLparm + Dx 
  }
  locsXaug <- sXaug[,-(parmcol_X),drop=FALSE]
  locszAug <- as.matrix(szAug-sXaug[,parmcol_X]*beta[parmcol_X]) # w z - w X beta....
  beta[-(parmcol_X)] <- get_from_MME(locsXaug,szAug=locszAug) 
  return(list(beta=beta)) # levQ is presumably always dense
}
