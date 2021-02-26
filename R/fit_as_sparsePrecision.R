.do_damped_WLS_spprec <- function(
  sXaug, zInfo, # fixed descriptors of system to be solved
  old_Vscaled_beta,
  oldAPHLs,
  APHLs_args,
  damping,
  dampingfactor=2, ## no need to change the input value
  ypos,off,GLMMbool,etaFix,
  updateW_ranefS_subarglist,
  wranefblob, seq_n_u_h, 
  ZAL_scaling, # locally fixed, "resident"; only changed in return value
  processed, 
  Trace,
  phi_est, #H_global_scale, 
  n_u_h, 
  ZAL, # unscaled one 
  which_LevMar_step,
  update_sXaug_constant_arglist,
  # promise rather than argument:
  low_pot=NULL,
  v_infer_args=NULL, # not null for beta optimization with v_in_b optimization i.e. in *some* .do_damped_WLS_outer() call
  stylefn, # in-loop stylefn for damped_WLS
  stylefn_v_out= .spaMM.data$options$stylefns$v_out_last, 
  stylefn_v_in= .spaMM.data$options$stylefns$v_in_last, ##
  outer) {
  ZAL_scaling <- 1 ## TAG: scaling for spprec
  if (outer) {
    trace <- max(0L,Trace-1L)
    stylefn_v <- stylefn_v_out
  } else {
    trace <- max(0L,Trace-2L)
    stylefn_v <- stylefn_v_in
  }
  if (processed$p_v_obj=="p_v" && which_LevMar_step!="v") { 
    objname <- "p_v" 
  } else { 
    objname <- APHLs_args$which <- "hlik"
  }

  oldlik <- oldAPHLs[[objname]]
  initdamping <- damping
  gainratio_grad <- zInfo$gainratio_grad
  gone_thru_mini <- (damping==1e-7)
  first_it <- TRUE
  prev_gainratio <- -Inf
  if (Trace && ! is.null(v_infer_args)) {
    cat(stylefn("[")) # cat(which_LevMar_step) #=> a substep of V_IN_B: "strict_v|b" or "b_&_v_in_b"
  } 
  loc_pot_tol <- attr(low_pot,"pot_tol") ## if low_pot is not NULL, we get this info.
  #v_total_maxit_mean <- v_infer_args$maxit.mean
  GLGLLM_const_w <- attr(processed$models,"GLGLLM_const_w")
  pot4improv <- NULL
  while ( TRUE ) { ## loop on damping; each iteration produce blue + ul-greens + yellow
    # if (trace && ! is.null(v_infer_args)) {
    #   cat(c(damping)) # => what is shown after [.v_h iter...]V_h IRLS...; .dampingfactor=* damping=
    # }
    if (processed$HL[1L]==1L) { ## ML fit 
      Vscaled_beta <- old_Vscaled_beta
      LM_z <- zInfo["scaled_grad"] # still a list, but for clarity, emphasizes that only this element is needed
      ## maximize p_v wrt beta only
      if ( which_LevMar_step=="v_b") { ## note tests on summand too !!!!!!!
        LevMarblob <- get_from_MME(sXaug=sXaug, which="LevMar_step", LM_z=LM_z, damping=damping)
        Vscaled_beta <- list(
          v_h=Vscaled_beta$v_h + LevMarblob$dVscaled, ##
          beta_eta=Vscaled_beta$beta_eta + LevMarblob$dbeta_eta
        )  
      } else if ( which_LevMar_step=="v") { ## v_h estimation given beta 
        LevMarblob <- get_from_MME(sXaug=sXaug, which="LevMar_step_v_h", LM_z=LM_z, damping=damping)
        Vscaled_beta$v_h <- Vscaled_beta$v_h + LevMarblob$dVscaled
      } else if ( which_LevMar_step=="strict_v|b") {
        # =: the case where we only fit v_h for the input beta_eta:
        #    only v_h in Vscaled_beta will be changed, by v_infer_args step below, 
        LevMarblob <- v_infer_args$LevMarblob ## LevMarblob$d... will be overwritten below
        #                                        is.null(LevMarblob) occurs in initial damping=Inf call.
      } else if ( which_LevMar_step=="b_&_v_in_b") {     
        LevMarblob <- get_from_MME(sXaug=sXaug, which="LevMar_step", LM_z=LM_z, damping=damping) # template
        dbeta_LevMarblob <- get_from_MME(sXaug=sXaug, which="LevMar_step_beta", LM_z=LM_z, damping=damping)
        Vscaled_beta$beta_eta <- Vscaled_beta$beta_eta + dbeta_LevMarblob$dbeta_eta
        # v_h in Vscaled_beta will be changed by v_infer_args step below.
        v_infer_args$LevMarblob <- get_from_MME(sXaug=sXaug, which="LevMar_step_v_h", LM_z=LM_z, damping=damping)
      } else if ( which_LevMar_step=="b_from_v_b") { 
        LevMarblob <- get_from_MME(sXaug=sXaug, which="LevMar_step", LM_z=LM_z, damping=damping)
        Vscaled_beta$beta_eta <- Vscaled_beta$beta_eta + dbeta_LevMarblob$dbeta_eta
      } else if ( which_LevMar_step=="b") { ## probably not used; get_from_MME() includes protection for pforpv=0
        LevMarblob <- get_from_MME(sXaug=sXaug, which="LevMar_step_beta", LM_z=LM_z, damping=damping)
        Vscaled_beta$beta_eta <- Vscaled_beta$beta_eta + dbeta_LevMarblob$dbeta_eta
      } 
      if ( which_LevMar_step!="v" &&  ! is.null(v_infer_args)) {
        if (Trace) cat(stylefn_v("["))
        #v_infer_args$maxit.mean <- ceiling(v_total_maxit_mean/5)
        v_h_blob <- .wrap_v_h_IRLS(v_h=Vscaled_beta$v_h , 
                                   beta_eta=Vscaled_beta$beta_eta, seq_n_u_h, GLMMbool, wranefblob, 
                                   processed$reserve$constant_u_h_v_h_args, updateW_ranefS_subarglist=updateW_ranefS_subarglist, 
                                   v_infer_args, Trace, IRLS_fn=".solve_v_h_IRLS_spprec")
        if (v_h_blob$break_info$IRLS_breakcond=="maxit") { # problematic failure: we need to do something
          #v_infer_args$maxit.mean <- ceiling(v_total_maxit_mean*4/5)
          psi_M <- rep(attr(processed$rand.families,"unique.psi_M"),diff(processed$cum_n_u_h))
          v_h_blob <- .wrap_v_h_IRLS(v_h=psi_M, 
                                     beta_eta=Vscaled_beta$beta_eta, seq_n_u_h, GLMMbool, wranefblob, 
                                     processed$reserve$constant_u_h_v_h_args, updateW_ranefS_subarglist=updateW_ranefS_subarglist, 
                                     v_infer_args, Trace, IRLS_fn=".solve_v_h_IRLS_spprec")
        } 
        # each underline green is a damping _loop_ not a damping step
        if (Trace) cat(stylefn_v("]"))
        if (trace)  { ## prints, at the level of the outer damped_WLS, the results of the v_h IRLS
          with(v_h_blob,cat(stylefn_v(paste0("v_h IRLS returns max(|grad|): v=",.prettify_num(break_info$maxs_grad[1L]), # grad when v_h IRLS exits
                             " beta=",.prettify_num(break_info$maxs_grad[2L]),
                             " after ",innerj, # number of iterations of v_h IRLS
                             " iter"
                             ))))
          break_info <- v_h_blob$break_info
          break_info$maxs_grad <- NULL
          for (st in names(break_info)) {
            if (is.numeric(stinfo <- break_info[[st]])) {
              cat(stylefn_v(paste0(", ",st,"=",.prettify_num(stinfo))))
            } else cat(stylefn_v(paste0(", ",st,"=",stinfo)))
          }
          cat(stylefn(";"))
        } # else if (Trace) cat(stylefn_v(".")) ## underline blue ## But no reason since no matrix factorization
        Vscaled_beta$v_h <- v_h_blob$v_h
        LevMarblob$dVscaled_beta <- list(
          v_h=Vscaled_beta$v_h - old_Vscaled_beta$v_h, ##
          beta_eta=Vscaled_beta$beta_eta -old_Vscaled_beta$beta_eta
        ) 
      }
    } else { ## joint hlik maximization
      LevMarblob <- get_from_MME(sXaug=sXaug, which="LevMar_step", LM_z=zInfo["scaled_grad"], damping=damping)
      Vscaled_beta <- list(
        v_h=old_Vscaled_beta$v_h + LevMarblob$dVscaled,
        beta_eta=old_Vscaled_beta$beta_eta + LevMarblob$dbeta_eta
      )  
    }
    v_h <- Vscaled_beta$v_h #### * ZAL_scaling =1 here
    eta <- off + drop(processed$AUGI0_ZX$X.pv %*% Vscaled_beta$beta_eta) + drop(ZAL %id*% v_h) ## length nobs 
    
    newmuetablob <- .muetafn(eta=eta,BinomialDen=processed$BinomialDen,processed=processed) 
    neww.resid <- .calc_w_resid(newmuetablob$GLMweights,phi_est)
    #newweight_X <- .calc_weight_X(neww.resid, H_global_scale) ## sqrt(s^2 W.resid)
    
    if (is.null(etaFix$v_h)) { 
      if (GLMMbool) {
        u_h <- v_h
        newwranefblob <- wranefblob ## keep input wranefblob since GLMM and lambda_est not changed
      } else {
        u_h_v_h_from_v_h_args <- c(processed$reserve$constant_u_h_v_h_args,list(v_h=v_h))
        u_h <- do.call(".u_h_v_h_from_v_h",u_h_v_h_from_v_h_args)
        if ( ! is.null(attr(u_h,"v_h"))) { ## second test = if constant_u_h_v_h_args$upper.v_h or $lower.v_h non NULL
          v_h <- attr(u_h,"v_h")
        }
        ## update functions u_h,v_h
        newwranefblob <- do.call(".updateW_ranefS",c(updateW_ranefS_subarglist, list(u_h=u_h,v_h=v_h)))
      } 
    } else newwranefblob <- wranefblob
    sXaug_arglist <- c(update_sXaug_constant_arglist,
                          list(w.ranef=newwranefblob$w.ranef, 
                               #weight_X=newweight_X,
                               w.resid=neww.resid))
    if ( ! GLMMbool) {newZAL_scaling <- 1}  ## TAG: scaling for spprec
    ####
    APHLs_args$dvdu <- newwranefblob$dvdu
    APHLs_args$u_h <- u_h 
    APHLs_args$mu <- newmuetablob$mu
    #
    if (processed$p_v_obj=="p_v" && which_LevMar_step!="v") { ## new damping -> new weights -> new expensive computation to evaluate p_v
      if (GLGLLM_const_w) {
        newsXaug <- NULL
        APHLs_args$sXaug <- sXaug
      } else {
        if (Trace) cat(stylefn(".")) # yellow in V_IN_B case
        newsXaug <- do.call(processed$AUGI0_ZX$envir$method, # ie, def_AUGI0_ZX_sparsePrecision
                            sXaug_arglist)
        APHLs_args$sXaug <- newsXaug
      }
    } else { ## sXaug_arglist will still be used after the loop !!!!!!!!!!!!!!!!!!!!
      APHLs_args$sXaug <- newsXaug <- NULL # to compute hlik, no expensive matrix computation.
    }
    newAPHLs <- do.call(".calc_APHLs_from_ZX", APHLs_args)
    #print(c(unlist(newAPHLs)))
    newlik <- unlist(newAPHLs[objname]) # keep name
    #
    if (damping==0L) {
      breakcond <- "damping=0"
      break
    }
    if ( which_LevMar_step=="strict_v|b") { ## whether rescue or not...
      breakcond <- "v|b_no_loop"
      break
    } #  =: single call to .calc_APHLs_from_ZX to only fit v_h for the input beta_eta.
    if (is.null(pot4improv)) { 
      switch(which_LevMar_step,
             "v_b"= {
               pot4improv <- sum(zInfo$m_grad_obj*unlist(get_from_MME(sXaug,szAug=zInfo["m_grad_obj"])))
               loc_pot_tol <- processed$spaMM_tol$b_pot_tol
             },
             "v"= {
               pot4improv <- get_from_MME(sXaug=sXaug, which="Mg_invH_g", B=gainratio_grad[seq_n_u_h])
               loc_pot_tol <- processed$spaMM_tol$v_pot_tol
             },
             "b"= { 
               B=zInfo$gainratio_grad[-seq_n_u_h]
               if (length(B)) {
                 pot4improv <- get_from_MME(sXaug=sXaug, which="Mg_invXtWX_g", B=B)
               } else pot4improv <- 0
               loc_pot_tol <- processed$spaMM_tol$b_pot_tol
             }, 
             "b_&_v_in_b"= {
               pot4improv <- sum(zInfo$m_grad_obj*unlist(get_from_MME(sXaug,szAug=zInfo["m_grad_obj"]))) 
               loc_pot_tol <- processed$spaMM_tol$b_pot_tol
             },
             "b_from_v_b"= {
               pot4improv <- get_from_MME(sXaug=sXaug, which="Mg_solve_g", B=zInfo$gainratio_grad) # fictitious, not implemented
               loc_pot_tol <- processed$spaMM_tol$b_pot_tol
             }
      )
      if (is.null(low_pot)) low_pot <- (pot4improv < loc_pot_tol) 
    } 
    if (low_pot) { # keeping the low_pot condition may be important for the "605" tests. We may always suppress it by the spaMM.options.
      very_low_pot <- (pot4improv < loc_pot_tol/10)
      if (processed$HL[1L]==1L && which_LevMar_step=="v_b") {
        v_pot4improv <- get_from_MME(sXaug=sXaug, which="Mg_invH_g", B=gainratio_grad[seq_n_u_h])
        breakcond <- structure("low_pot", pot4improv=pot4improv, very_low_pot=very_low_pot,
                               no_overfit = (v_pot4improv < processed$spaMM_tol$v_pot_tol))
      } else breakcond <- structure("low_pot", pot4improv=pot4improv, very_low_pot=very_low_pot)
      break
    } 
    ## ELSE
    gainratio <- (newlik!=-Inf) ## -Inf occurred in binary probit with extreme eta... 
    if (gainratio) { 
      summand <- .calc_summand_gainratio_spprec(processed, which_LevMar_step, LevMarblob, seq_n_u_h, ZAL_scaling=1, gainratio_grad)
      ## The two terms of the summand should be positive. In part. dv_h_beta*gainratio_grad should be positive. 
      ## However, numerical error may lead to <0 or even -Inf
      ## Further, if there are both -Inf and +Inf elements the sum is NaN.
      summand[summand<0] <- 0
      denomGainratio <- sum(summand)
      dlogL <- newlik-oldlik
      conv_logL <- abs(dlogL)/(1+abs(newlik))
      gainratio <- 2*dlogL/denomGainratio ## cf computation in MadsenNT04 below 3.14, but generalized for D' != I ; numerator is reversed for maximization
      # but if gradient is practically zero and damping  ~0 we may not wish to compare ~0 to ~0...
      #  which is why we break, rather than stop, on (damping>1e100). MadsenNT04 have other stopping crits 
    } else { ## 2017/10/16 a patch to prevent a stop in this case, but covers up dubious computations (FIXME)
      newlik <- -.Machine$double.xmax
      dlogL <- newlik-oldlik
      conv_logL <- abs(dlogL)/(1+abs(newlik))
      denomGainratio <- Inf # bc denomGainratio can be tested below 
    }
    if (trace) {
      cat(stylefn(paste(" dampingfactor=",dampingfactor,#"innerj=",innerj,
                           "damping=",damping,"gainratio=",gainratio,"oldlik=",oldlik,"newlik=",newlik)))
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
      if (div_gainratio < -max(LevMarblob$dampDpD)/(1e05*damping)) { 
        breakcond <- "div_gain" # used to switch from "V_IN_B" to "strict_v|b"
        break 
      } 
    }
    if (gainratio > 0) { # gainratio may be always negative if initial ranefs better optimize logL than the correct solution does.
      ## cf Madsen-Nielsen-Tingleff again, and as in levmar library by Lourakis
      damping <- damping * max(1/3,1-(2*gainratio-1)^3) # gainratio->0 factor-> 2; gainratio->1 factor->1/3
      breakcond <- "OK_gain"
      break 
    }
    if ( dampingfactor>4 ## ie at least 2 iteration of the while() => prev_conv_logL is available
                && conv_logL <1e-8 && abs(prev_conv_logL) <1e-8 
                && gone_thru_mini ## has gone through 1e-7, including the cases where we started with too high damping 
    ) { 
      breakcond <- "stuck_obj"
      break   ##   cases were we do not expect any significant improvement
    } 
    if ( (! gone_thru_mini) &&  # condition to avoid infinite loop when gainratio remains <0 for arbitrarily large damping
         #conv_logL < 1e-9 # not useful as this may mean e were borderline making progress
         # in which case we should damp further! low_pot is not instruction bc it must be FALSE when we reach here. Another indicator of potential is:
         denomGainratio<loc_pot_tol/10 # lax condition here leads to poor perf.
    ) { 
      damping <- 1e-7
      dampingfactor <- 2
      gone_thru_mini <- TRUE
      prev_gainratio <- -Inf
      if (trace) cat("-")
      # and continue 
    } else { ## other UNsuccessful step
      prev_gainratio <- gainratio
      prev_conv_logL <- conv_logL
      damping <- damping*dampingfactor
      dampingfactor <- dampingfactor*2
      if (damping>1e100) { # endpoint for large negative gainratio (i.e. v overfit as starting point)
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
  if (is.null(newsXaug)) { ## which means that hlik is the local objective.
    # For HL11, p_v will be used as oldAPHLs in the next call to .do_damped_WLS_outer() in an alternating algo;
    #   and sXaug may be needed to compute sscaled in .solve_v_h_IRLS()
    # For PQL fits newsXaug has not been needed in the damping loop but will be needed after exiting this fn
    #   (e.g., for its next call -> LevMarblob <- get_from_MME(sXaug=sXaug, which="LevMar_step", LMrhs=zInfo$scaled_grad, damping=damping))
    if (GLGLLM_const_w) {
      APHLs_args$sXaug <- newsXaug <- sXaug
    } else {
      if (Trace) { 
        if (processed$p_v_obj=="p_v") { # v estimation within HL11
          cat(stylefn_v("."))
        } else  cat(stylefn(".")) # PQL/L, vb extimation
      }
      newsXaug <- do.call(processed$AUGI0_ZX$envir$method, sXaug_arglist)
      APHLs_args$sXaug <- newsXaug
    } 
    APHLs_args$which <- processed$p_v_obj # "p_v" # 
    newAPHLs <- do.call(".calc_APHLs_from_ZX", APHLs_args) 
  }
  RESU <- list(lik=newlik, APHLs=newAPHLs, damping=damping, sXaug=newsXaug,
               # fitted=fitted, ## FIXME: removed so that no shortcut a la Bates in calc_APHLs_from_ZX; reimplement the shorcut in that fn?
               eta=newmuetablob$sane_eta, muetablob=newmuetablob, wranefblob=newwranefblob,
               breakcond=breakcond,
               v_h=v_h, u_h=u_h, w.resid=neww.resid) # newweight_X does not exists and not needed in spprec
  if ( ! first_it) { # if not break in first iteration
    RESU$conv_logL_not_first_it <- conv_logL
  }
  if ( ! GLMMbool ) {
    RESU$ZAL_scaling <- newZAL_scaling
    # RESU$Xscal <- newXscal ## does not exist and presumably not needed.
    Vscaled_beta$v_h <- v_h/newZAL_scaling ## represent solution in new scaling...
  } 
  RESU$Vscaled_beta <- Vscaled_beta 
  return(RESU)
}

#copies to allow independent debug()ing
.do_damped_WLS_v_in_b_spprec <- .do_damped_WLS_spprec
.do_damped_WLS_outer_spprec <- .do_damped_WLS_spprec

# processed, Trace


.WLS_substitute_spprec <- function(update_sXaug_constant_arglist, Vscaled_beta, off, etaFix, mod_attr, 
                                   updateW_ranefS_subarglist, ZAL, 
                            processed, phi_est,
                            wranefblob, Trace,stylefn) {
  
  # Vscaled_beta must have been provided by somethin else than damped_WLS_blob
  # drop, not as.vector(): names are then those of (final) eta and mu -> used by predict() when no new data
  eta <- off + drop(processed$AUGI0_ZX$X.pv %*% Vscaled_beta$beta_eta) + drop(ZAL %id*% Vscaled_beta$v_h)
  RESU <- list()
  if (is.null(etaFix$v_h)) { 
    v_h <- Vscaled_beta$v_h ## * ZAL_scaling (=1)
    if ( mod_attr$GLMMbool ) {
      RESU$u_h <- RESU$v_h <- v_h ## keep input wranefblob since lambda_est not changed
    } else {
      u_h_v_h_from_v_h_args <- c(processed$reserve$constant_u_h_v_h_args,list(v_h=v_h))
      RESU$u_h <- u_h <- do.call(".u_h_v_h_from_v_h",u_h_v_h_from_v_h_args)
      if ( ! is.null(attr(u_h,"v_h"))) { ## second test = if constant_u_h_v_h_args$upper.v_h or $lower.v_h non NULL
        v_h <- attr(u_h,"v_h")
      }
      RESU$v_h <- v_h
      ## update functions u_h,v_h
      RESU$wranefblob <- wranefblob <- do.call(".updateW_ranefS",c(updateW_ranefS_subarglist, list(u_h=u_h,v_h=v_h)))
      #if ( ! mod_attr$GLMMbool) { 
        RESU$ZAL_scaling <- 1 ## TAG: scaling for spprec
      # } 
    }
  }
  RESU$muetablob <- muetablob <- .muetafn(eta=eta,BinomialDen=processed$BinomialDen,processed=processed) 
  if ( (! mod_attr$LLM_const_w) && (! mod_attr$GLGLLM_const_w) ) {
    RESU$w.resid <- .calc_w_resid(muetablob$GLMweights,phi_est)
    sXaug_arglist <- c(update_sXaug_constant_arglist, # contained H_global_scale but not longer so
                           list(w.ranef=wranefblob$w.ranef, 
                                #weight_X=weight_X, 
                                w.resid=RESU$w.resid))
    if (Trace) cat(stylefn("."))
    RESU$sXaug <- do.call(processed$AUGI0_ZX$envir$method, # ie, def_AUGI0_ZX_sparsePrecision
                          sXaug_arglist)
  } ## ergo sXaug is not updated for LMM (no need to)
  return(RESU) ## contains only updated quantities
}

.solve_IRLS_as_spprec <- 
  function(
           ZAL, y, 
           n_u_h, 
           #H_global_scale, 
           lambda_est, muetablob=NULL, off, maxit.mean, etaFix,
           wranefblob, processed,
           ## for ! LMM
           phi_est, 
           ## supplement for LevenbergM
           beta_eta,
           ## supplement for ! GLMM
           u_h, v_h, w.resid=NULL, 
           for_init_z_args,
           for_intervals,
           ##
           corrPars, # corrPars needed together with adjMatrix to define Qmat
           verbose,
           LevM_HL11_method=.spaMM.data$options$LevM_HL11_method
  ) {
    trace <- verbose["TRACE"]
    if (trace) {
      cat(">") 
      if (verbose["trace"]) cat(.pretty_summ_lambda(lambda_est,processed))
    }
  pforpv <- ncol(processed$AUGI0_ZX$X.pv)
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
    which_LevMar_step <- default_b_step <- LevM_HL11_method[["b_step"]] 
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
    constant_zAug_args <- list(n_u_h=n_u_h, nobs=nobs, pforpv=pforpv, y=y, off=off, ZAL=ZAL, processed=processed)
    if ( ! GLMMbool) {
      constant_init_z_args <- c(list(lcrandfamfam=lcrandfamfam, nobs=nobs, lambda_est=lambda_est, ZAL=ZAL),  
                                # fit_as_ZX args specific for ! GLMM:
                                for_init_z_args,
                                #
                                mget(c("cum_n_u_h","rand.families"),envir=processed))
      updateW_ranefS_subarglist <- processed$reserve$W_ranefS_constant_args
      updateW_ranefS_subarglist$lambda <- lambda_est
    } 
  } 
  
  ##### initial sXaug
  ZAL_scaling <- 1  ## TAG: scaling for spprec
  if (is.null(muetablob)) { ## NULL input eta allows NULL input muetablob
    eta <- off + drop(processed$AUGI0_ZX$X.pv %*% beta_eta) + drop(ZAL %id*% v_h) 
    muetablob <- .muetafn(eta=eta,BinomialDen=processed$BinomialDen,processed=processed) 
  }
  ## varies within loop if ! LMM since at least the GLMweights in w.resid change
  if ( is.null(w.resid) ) w.resid <- .calc_w_resid(muetablob$GLMweights,phi_est)
  ## needs adjMatrix and corrPars to define Qmat
  update_sXaug_constant_arglist <- list(AUGI0_ZX=processed$AUGI0_ZX, corrPars=corrPars, 
                                        cum_n_u_h=processed$cum_n_u_h #,H_global_scale=H_global_scale
                                        ) 
  #weight_X <- .calc_weight_X(w.resid, H_global_scale) ## sqrt(s^2 W.resid)
  if (trace) {
    stylefn <- switch(which_LevMar_step,
                      v=.spaMM.data$options$stylefns$vloop,
                      V_IN_B=.spaMM.data$options$stylefns$v_in_loop,
                      .spaMM.data$options$stylefns$betaloop )
    if (LevenbergM) cat("LM")
    cat(stylefn("."))
  }
  sXaug <- do.call(processed$AUGI0_ZX$envir$method, # ie, def_AUGI0_ZX_sparsePrecision
                   c(update_sXaug_constant_arglist,
                     list(w.ranef=wranefblob$w.ranef, 
                          #weight_X=weight_X, 
                          w.resid=w.resid)))
  if ( ! is.null(for_intervals)) {
    Vscaled_beta <- list(v_h=v_h/ZAL_scaling, beta_eta=for_intervals$beta_eta)
    fixefobjfn <- names(for_intervals$fixeflik)
  } else {
    Vscaled_beta <- list(v_h=v_h/ZAL_scaling,beta_eta=beta_eta)
  } 
  # to be evaluated once when it becomes needed:
  delayedAssign("constant_v_infer_args", list( # ultimately for the .solve_v_h_IRLS_spprec() call
    X.pv=processed$AUGI0_ZX$X.pv, 
    ZAL=ZAL, y=y, n_u_h=n_u_h, #H_global_scale=H_global_scale,
    lambda_est=lambda_est, off=off,maxit.mean=maxit.mean,etaFix=etaFix,
    processed=processed, phi_est=phi_est, for_init_z_args=for_init_z_args,
    trace=trace, corrPars=corrPars, dampings_env=dampings_env))
  ## Loop controls:
  allow_LM_restart <- ( ! LMMbool && ! LevenbergM && is.null(for_intervals) && is.na(processed$LevenbergM["user_LM"]) )
  if (allow_LM_restart) {
    keep_init <- new.env()
    names_keep <- c("sXaug","wranefblob","muetablob","u_h","w.resid","v_h","beta_eta","Vscaled_beta")
    for (st in names_keep) keep_init[[st]] <- environment()[[st]]
  }
  LMcond <- - 10. # also for hlik LM algo
  ################ L O O P ##############
  for (innerj in 1:maxit.mean) {
    if( ! LevenbergM && allow_LM_restart) { ## FIXME the next step improvement would be 
      #  ./. to keep track of lowest lambda that created problem and use LM by default then
      if (innerj>3) {
        LMcond <- LMcond + mean(abs_d_relV_beta/(old_abs_d_relV_beta+1e-8))
        ## cat(mean(abs_d_relV_beta/old_abs_d_relV_beta)," ")
        # cat(LMcond/innerj," ")
        if (LMcond/innerj>0.5) {
          if (trace) cat("!LM") # ie, LevenbergM!
          for (st in names_keep) assign(st,keep_init[[st]])
          LevenbergM <- TRUE
          # Vscaled_beta included in keep_init hence assigned by assign(st,keep_init[[st]])
          dampings_env <- list2env(.spaMM.data$options$spaMM_tol$dampings_env_v)
          damped_WLS_blob <- NULL
          allow_LM_restart <- FALSE 
          if ( which_LevMar_step=="v_b") { 
            ## The LevM.negbin test finds "strict_v|b" poorer than "V_IN_B"" (note some divergent p_v's)  which led to:
            # which_LevMar_step <- "V_IN_B" ## not modified by if (... ! is.null(damped_WLS_blob) ...) before being used.
            # but the  optim_LevM's (update(br$fullfit,fixed=... test shows one should keep using "v_b" here.
            # otherwise "!LM" differs (and is poorer) from LevM=TRUE (which indeed starts from "v_b")
            # However, it's not clear why "V_IN_B" is poorer (and it's not the step on which the loop terminates)  => (there is a fixme on the fit_as_ZX version...) 
            from_good_v_b <- FALSE
            rescue_thr <- processed$spaMM_tol$rescue_thr
            rescue_nbr <- 0L
            prev1_rescued <- FALSE
          }
        }
      }
      if (innerj>2) old_abs_d_relV_beta <- abs_d_relV_beta
    }
    ##### get the lik of the current state
    if ( ! is.null(for_intervals)) {
      loc_logLik_args <- list(sXaug=sXaug, processed=processed, phi_est=phi_est,
                              lambda_est=lambda_est, dvdu=wranefblob$dvdu, u_h=u_h, mu=muetablob$mu)
      oldlik <- unlist(do.call(".calc_APHLs_from_ZX",loc_logLik_args)[fixefobjfn]) # unlist keeps name
    } else if (LevenbergM) { ## then logL is necessary to check for increase
      if (is.null(damped_WLS_blob)) {
        oldAPHLs <- .calc_APHLs_from_ZX(sXaug=sXaug, processed=processed, phi_est=phi_est, which=processed$p_v_obj, 
                                        lambda_est=lambda_est, dvdu=wranefblob$dvdu, u_h=u_h, mu=muetablob$mu)
      } else { ## Levenberg and innerj>1
        oldAPHLs <- damped_WLS_blob$APHLs
      }
    } 
    #####
    
    ##### RHS
    if (LMMbool) {
      zInfo <- list(z2=NULL,z1=y-off,sscaled=0)
      zInfo$z1_eta <- z1_sscaled_eta <- y-muetablob$sane_eta
    } else {
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
                               init_z_args=init_z_args) )
      zInfo <- do.call(".calc_zAug_not_LMM",calc_zAug_args)
      if (GLMMbool) zInfo$z2 <- NULL
      etamo <- muetablob$sane_eta - off
      zInfo$z1_eta <- zInfo$z1- etamo 
      z1_sscaled_eta <- zInfo$z1_sscaled - etamo # augz[-seq_n_u_h]- etamo # z_1-sscaled-etamo
    }
    ## keep name 'w'zAug to emphasize the distinct weightings  of zaug and Xaug (should have been so everywhere)
    #####
    ##### improved  Vscaled_beta   
    ## he solver uses a 'beta first" approach to solveing for d_beta and d_v... even for LMM 
    if (GLMMbool) {
      zInfo$dlogfvdv <-  - v_h * wranefblob$w.ranef
    } else zInfo$dlogfvdv <- (zInfo$z2 - v_h) * wranefblob$w.ranef
    ## $gainratio_grad is the rhs in the direct solution of the full system (by chol2inv in the dense QR case)
    ## and thus propto the gradient of the objective.
    ## The 'betaFirst' algo in sparsePrecision does not use it, only the gainratio code does  
    ## the gradient for -p_v (or -h), independent of the scaling
    zInfo$m_grad_obj <- c( ## drop() avoids c(Matrix..); attr(sXaug,"w.resid") correct for truncated models too.
      m_grad_v <- drop(.crossprod(ZAL, attr(sXaug,"w.resid") * zInfo$z1_eta) +zInfo$dlogfvdv), # Z'W(z_1-eta) + dlogfvfv
      drop(.crossprod(processed$AUGI0_ZX$X.pv, attr(sXaug,"w.resid") * z1_sscaled_eta)) # X'W(z_1-sscaled-eta) 
    )
    if (LevenbergM) {
      zInfo$gainratio_grad <- zInfo$m_grad_obj 
      zInfo$scaled_grad <- zInfo$m_grad_obj
    }
    #
    if ( ! is.null(for_intervals)) {
      currentDy <- (for_intervals$fixeflik-oldlik)
      if (LMMbool) etamo <-muetablob$sane_eta - off
      intervalBlob <- .intervalStep_spprec(old_v_h_beta=Vscaled_beta,
                                       sXaug=sXaug,zInfo=zInfo,
                                       for_intervals=for_intervals,
                                       currentlik=oldlik,currentDy=currentDy,
                                       ZAL=ZAL, etamo=etamo)
      damped_WLS_blob <- NULL
      Vscaled_beta <- intervalBlob$v_h_beta
    } else if (LevenbergM) {
      if (trace>1L) {
        stylefn <- switch(which_LevMar_step,
                       v=.spaMM.data$options$stylefns$vloop,
                       V_IN_B=.spaMM.data$options$stylefns$v_in_loop,
                       .spaMM.data$options$stylefns$betaloop )
        maxs_grad <- c(max(abs(zInfo$m_grad_obj[seq_n_u_h])),max(abs(zInfo$m_grad_obj[-seq_n_u_h])))
        cat(stylefn("iter=",innerj,", max(|grad|): v=",maxs_grad[1L],"beta=",maxs_grad[2L],";"))
      }
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
          if (prev1_rescued || rescue_nbr > rescue_thr["V_IN_B"]) {
            which_LevMar_step <- "V_IN_B"
          } else which_LevMar_step <- "v_b" 
        } else if (which_LevMar_step=="v_b" || which_LevMar_step=="b_from_v_b" ) { 
          if (rescue_nbr > rescue_thr["strictv"]  &&  #rescue has been previously needed in the outer loop
              damped_WLS_blob$breakcond != "low_pot" ) { # i.e. if OK_gain (other cases would lead to the previous "rescue" case)
            which_LevMar_step <- "strict_v|b" # yet we play safer if we know a problem occurred previously
            from_good_v_b <- TRUE
          # } else if (max(abs(m_grad_v)) > max(abs(old_m_grad_v))) which_LevMar_step <- "v" # test is fausse bonne idee...
          } else which_LevMar_step <- "v"
        } else if (which_LevMar_step=="v") {
          if (damped_WLS_blob$breakcond == "low_pot") { ## LevMar apparently maximized h wrt v after several iterations
            #cat(damped_WLS_blob$breakcond)
            old_relV_beta <- relV_beta ## serves to assess convergence !!! which is thus dependent on condition ( hlik_stuck || ! need_v_step)
            which_LevMar_step <- default_b_step # We should not reach this line when RHS is "v_in_b"
          } else {
            ## v_h estimation not yet converged, continue with it 
            # But this means that "stuck_obj" on a "v" step must have dealt with elsewhere, otherwise we are looping on stuck_obj
            # => special rescue case.
          }
        } else if (which_LevMar_step=="strict_v|b") {
          old_relV_beta <- relV_beta 
          if (from_good_v_b) { # strictv was called after a v_b (which implies that rescue was not called)
            which_LevMar_step <- "v_b"
          } else { # strictv was called after V_IN_B (!)
            if (rescue_nbr > rescue_thr["re_V_IN_B"]) {
              which_LevMar_step <- "V_IN_B"
            } else which_LevMar_step <- "v_b"  
          }
          #which_LevMar_step <- "V_IN_B" # sequence "v_b" -> "low_pot" -> "strict_v|b"
        } else if (which_LevMar_step=="V_IN_B") { 
          breakcond <- damped_WLS_blob$breakcond
          if (breakcond=="stuck_obj" || breakcond=="div_gain") {
            which_LevMar_step <- "strict_v|b" # the call to "strict_v|b" may seem odd but results in clean optim
            #If we did that in the wrap... then we would next compare two identical "strict_v|b" 
            from_good_v_b <- FALSE
          } else {
            old_relV_beta <- relV_beta 
            if (rescue_nbr > rescue_thr["re_V_IN_B"]) {
              which_LevMar_step <- "V_IN_B"
            } else which_LevMar_step <- "v_b" 
          }
        } else if (default_b_step=="v_in_b") { # presumably not used
          old_relV_beta <- relV_beta 
        } else { ## "b" or any unanticipated case # presumably not used 
          # as v_b cas:
          if (rescue_nbr > rescue_thr["strictv"] &&  #rescue has been previously needed in the outer loop
              damped_WLS_blob$breakcond != "low_pot") { #rescue has been previously needed in the outer loop
            which_LevMar_step <- "strict_v|b" # yet we play safer if we know a problem ocurred previously
            from_good_v_b <- TRUE
          } else which_LevMar_step <- "v" # yet we play safer if we know a problem ocurred previously
        }
        prev1_rescued <- just_rescued 
      } 
      new_damping <- .get_new_damping(dampings_env$v[[which_LevMar_step]], which_LevMar_step)
      damped_WLS_blob <- .wrap_do_damped_WLS_outer(
        damped_WLS_fn = .do_damped_WLS_outer_spprec,
        LevM_HL11_method=LevM_HL11_method, # contains the rescue_thr options => any possibility to simplify arguments ?
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
        updateW_ranefS_subarglist=updateW_ranefS_subarglist,
        wranefblob=wranefblob,seq_n_u_h=seq_n_u_h,
        update_sXaug_constant_arglist=update_sXaug_constant_arglist,
        #H_global_scale=H_global_scale,     
        ZAL_scaling= ZAL_scaling, 
        processed=processed, 
        phi_est=phi_est, n_u_h=n_u_h, ZAL=ZAL,
        constant_v_infer_args=constant_v_infer_args,
        looseness= if ( is.null(damped_WLS_blob) ||  ## start strict
                        new_damping>1e-7) {## use strict when there are trace of difficulties (in particular, failure to improve) 
          1 } else {processed$spaMM_tol$loose_fac},
        low_pot=NULL ## explicit for clarity, but its the default
      ) 
      #old_m_grad_v <- m_grad_v
      dampings_env$v[[attr(damped_WLS_blob,"step")]] <- damped_WLS_blob$damping
      ## LevM PQL
      if (! is_HL1_1) {
        if (damped_WLS_blob$lik < oldAPHLs$hlik) { ## if LevM step failed to find a damping that increases the hlik
          # Tis should occur only bc of (1) numerically challenging conditions e.g mu close to bounds; or (2) optimum has been 
          # found and floating point innacurracies matter.  We try to exclude the second case by the following test: 
          if ( ! ((breakcond <- damped_WLS_blob$breakcond)=="low_pot" && attr(breakcond,"very_low_pot"))) {
            damped_WLS_blob <- NULL
            beta_cov_info <- .calc_beta_cov_info_spprec(X.pv = sXaug$AUGI0_ZX$X.pv,envir = sXaug$BLOB,w.resid = w.resid)
            if ((current_kappa <- kappa(tcrossprod(beta_cov_info$tcrossfac_beta_v_cov)))>1e08) {
              processed$envir$PQLdivinfo$high_kappa$ranFixes <- c(processed$envir$PQLdivinfo$high_kappa$ranFixes,
                                                                  list(processed$envir$ranFix))
            } else processed$envir$PQLdivinfo$unknown$ranFixes  <- c(processed$envir$PQLdivinfo$unknown$ranFixes,
                                                                     list(processed$envir$ranFix))
            dVscaled_beta <- get_from_MME(sXaug,szAug=zInfo) ################### FIT
            Vscaled_beta <- list(v_h=Vscaled_beta$v_h+dVscaled_beta$dv_h,
                                 beta_eta=Vscaled_beta$beta_eta+dVscaled_beta$dbeta_eta)
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
    } else { ## IRLS: always accept new v_h_beta
      dVscaled_beta <- get_from_MME(sXaug,szAug=zInfo) ################### FIT
      Vscaled_beta <- list(v_h=Vscaled_beta$v_h+dVscaled_beta$dv_h,
                           beta_eta=Vscaled_beta$beta_eta+dVscaled_beta$dbeta_eta)
      damped_WLS_blob <- NULL
    }
    if (trace>5L) .prompt()
    ##### Everything that is needed for 
    #  (1) assessment of convergence: c(v_h*sqrt(wranefblob$w.ranef),beta_eta)
    #  (2) all return elements are updated as function of the latest Vscaled_beta 
    #                               (itself possible updated to new scaling by the following assign()'s).
    #      In particular We need muetablob and (if ! LMM) sXaug, hence a lot of stuff.
    #  Hence, the following code is useful whether a break occurs or not. 
    if ( is.null(damped_WLS_blob) ) { ## fits nothing, but updates variables in case of standard IRLS, or of intervals
      WLS_blob <- .WLS_substitute_spprec(update_sXaug_constant_arglist, Vscaled_beta, off, etaFix, mod_attr=mod_attr, 
                                         updateW_ranefS_subarglist=updateW_ranefS_subarglist, ZAL, 
                                         processed, phi_est,
                                         wranefblob, Trace=trace, stylefn)
      for (st in names(WLS_blob)) assign(st,WLS_blob[[st]]) 
      if ( ! LMMbool && fpot_cond && is.null(for_intervals)) {
        Mg_solve_g <- sum(.unlist(dVscaled_beta)*zInfo$m_grad_obj) # using old m_grad_obj
        if (Mg_solve_g < fpot_tol) break
      } 
    } else {
      for (st in intersect(names(damped_WLS_blob), # sXaug (and the weights) need not be present if(GLGLLM_const_w) 
                           c("w.resid", ## !important! cf test-adjacency-corrMatrix.R
                             "sXaug","Vscaled_beta","wranefblob","v_h","u_h","muetablob"))) assign(st,damped_WLS_blob[[st]])
      if ( ! GLMMbool ) ZAL_scaling <- damped_WLS_blob$ZAL_scaling ## cf 'TAG:', this line in case I decide to use a ZAL_scaling !=1 later
    }
    #  At this point all return elements are updated as function of the latest Vscaled_beta.
    #  In particular We need muetablob and (if ! LMM) sXaug, hence a lot of stuff.
    #####
    beta_eta <- Vscaled_beta$beta_eta
    ##### assessment of convergence
    if (innerj==maxit.mean) {
      if (maxit.mean>1L) {
        if (LevenbergM) processed$LevenbergM["LM_start"] <- TRUE
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
      relV_beta <- c(v_h*sqrt(wranefblob$w.ranef),beta_eta)  ## convergence on v_h relative to sqrt(lambda), more exactly for Gaussian
      abs_d_relV_beta <- abs(relV_beta - old_relV_beta) ## for ML, comparison between estimates when ( hlik_stuck || ! need_v_step )
      if (is_HL1_1 && LevenbergM) {
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
        break
      }
      # More ad hoc breaks for cases where the coefficients keep moving although the total potential is low:
      if ( ! is.null(damped_WLS_blob) ) {
        if (is_HL1_1) {
          if ( which_LevMar_step=="v_b" && damped_WLS_blob$breakcond=="low_pot" && attr(damped_WLS_blob$breakcond,"no_overfit")) {
            break # motivated by BINARYboot (e.g., replicate 362)
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
      processed$envir$confint_best$beta_eta <- .unscale(sXaug$AUGI0_ZX$X.pv, Vscaled_beta[n_u_h+seq_len(pforpv)])
      processed$envir$confint_best$ranPars <- for_intervals$ranFix
      processed$envir$confint_best$ranPars$lambda <- unique(lambda_est) # ugly but get_inits_from_fit() cannot be used in this context...
    }
  }
  if (trace>1L && (LevenbergM))  {
    stylefn <- .spaMM.data$options$stylefns$betalast
    maxs_grad <- c(max(abs(zInfo$m_grad_obj[seq_n_u_h])),max(abs(zInfo$m_grad_obj[-seq_n_u_h])))
    cat(stylefn("iter=",innerj,", max(|grad|): v=",maxs_grad[1L],"beta=",maxs_grad[2L],";"))
  }
  names(beta_eta) <- colnames(processed$AUGI0_ZX$X.pv)
  RESU <- list(sXaug=sXaug, 
               ## used by calc_APHLs_from_ZX
               #fitted=fitted, ## FIXME: removed so that no shortcut a la Bates in calc_APHLs_from_ZX; reimplement the shorcut in that fn?  
               #weight_X=NA,
               nobs=nobs, pforpv=pforpv, seq_n_u_h=seq_n_u_h, u_h=u_h, 
               muetablob=muetablob, 
               lambda_est=lambda_est,
               phi_est=phi_est,
               ## used by other code
               beta_eta=beta_eta, w.resid=w.resid, wranefblob=wranefblob, 
               v_h=v_h, eta=muetablob$sane_eta, innerj=innerj)
  ## for diagnostic purposes
  if ( ! LMMbool ) RESU$m_grad_obj <- zInfo$gainratio_grad ## frm ZInfo bc other copies of m_grad_obj may be missing if not LevenbergM
  return(RESU)
} 

.intervalStep_spprec <- function(old_v_h_beta,sXaug,zInfo,currentlik,for_intervals,currentDy,ZAL, etamo) { 
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
    targetDy <- (for_intervals$fixeflik-for_intervals$targetlik)
    Dx <- currentDx*sqrt(targetDy/currentDy)
    #cat(currentDx," ",targetDy," ",currentDy," ",Dx,"\n")
    ## pb is if Dx=0 , Dx'=0... and Dx=0 can occur while p_v is still far from the target, because other params have not converged.
    ## (fixme) patch:
    if (currentDy<targetDy) { ## we are close to the ML: we extrapolate a bit more confidently
      min_abs_Dx <- for_intervals$asympto_abs_Dparm/1000
    } else min_abs_Dx <- 1e-6 ## far from ML: more cautious move our of Dx=0
    Dx <- sign(currentDx)*max(abs(Dx),min_abs_Dx)
    v_h_beta_vec[parmcol_ZX] <- for_intervals$MLparm + Dx 
  }
  ## gradient changes: (the get_from_MME() call only uses m_grad_obj to get rhs)
  # scratch code that led o this is in devel_confint_spprec.R
  parmcol_X <- for_intervals$parmcol_X
  off_newparm <- sXaug$AUGI0_ZX$X.pv[,parmcol_X]*v_h_beta_vec[parmcol_ZX] # function of NEW tentative value of MLparm
  # so it does not give the difference between etamo and i_etamo, which depends on OLD value of MLparm
  oovb <- old_v_h_beta_vec[-(parmcol_ZX)]
  seq_n_u_h <- seq_along(old_v_h_beta$v_h)
  i_etamo <- drop(ZAL %*% oovb[seq_n_u_h] +
                   sXaug$AUGI0_ZX$X.pv[,-(parmcol_X),drop=FALSE] %*%  oovb[-seq_n_u_h]) # 
  #
  rhs <- attr(sXaug, "w.resid")*(etamo-off_newparm- i_etamo)
  zInfo$m_grad_obj <- zInfo$m_grad_obj[-parmcol_ZX]+
    c(drop(crossprod(ZAL,rhs)),drop(crossprod(sXaug$AUGI0_ZX$X.pv[,-parmcol_X,drop=FALSE],rhs))) # tjrs zInfoo$m_grad_obj
  ## sXaug changes: 
  int_AUGI0_ZX <- as.list(sXaug$AUGI0_ZX) ## as.list() local copy avoids global modifs of original sXaug$AUGI0_ZX envir
  int_AUGI0_ZX$X.pv <- int_AUGI0_ZX$X.pv[,-(parmcol_X),drop=FALSE]
  int_AUGI0_ZX$ZeroBlock <- int_AUGI0_ZX$ZeroBlock[,-(parmcol_X),drop=FALSE]
  # Assuming sXaug not being an envir, sXaug$AUGI0_ZX <- int_AUGI0_ZX would be correct, creating a local copy of sXaug, 
  #  but less explicit than the following: 
  int_sXaug <- sXaug ## Assuming sXaug not being an envir, then the following does not affect sXaug.
  int_sXaug$AUGI0_ZX <- int_AUGI0_ZX ## Replaces an envir by a list i nthe local copy; 
  # ?__F I X M E__? it's unsafe as $BLOB is left unchanged (and still belong to the original object) while it might contain info reused to solve the system:
  # Theis attr(.,"pforpv") is used to determined whether there are beta coeffs to compute in .AUGI0_ZX_sparsePrecision()
  # int_sXaug$BLOB <- list2env(list(), parent=environment(.AUGI0_ZX_sparsePrecision))
  attr(int_sXaug,"pforpv") <- attr(int_sXaug,"pforpv") - length(parmcol_X) # a priori  -1
  v_h_beta_vec[-(parmcol_ZX)] <- old_v_h_beta_vec[-(parmcol_ZX)]+ unlist(get_from_MME(int_sXaug,szAug=zInfo)) 
  return(list(v_h_beta=relist(v_h_beta_vec,old_v_h_beta))) 
}

