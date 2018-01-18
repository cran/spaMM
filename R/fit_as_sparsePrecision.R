.do_damped_WLS_spprec <- function(sXaug, zInfo, # fixed descriptors of system to be solved
                                 old_v_h_beta, old_eta,
                                 oldAPHLs,
                                 APHLs_args,
                                 damping,
                                 dampingfactor=2, ## no need to change the input value
                                 ypos,off,GLMMbool,etaFix,constant_u_h_v_h_args,
                                 updateW_ranefS_constant_arglist,
                                 wranefblob, update_sXaug_constant_arglist,
                                 seq_n_u_h,
                                 processed, trace=FALSE,
                                 phi_est, n_u_h, ZAL,
                                 which_LevMar_step
) {
  initdamping <- damping
  gainratio_grad <- zInfo$m_grad_obj
  while ( TRUE ) {
    if (processed$HL[1L]==1L) { ## ML fit 
      if ( which_LevMar_step=="v_b") { ## note testes on summand too !!!!!!!
        LevMarblob <- get_from_MME(sXaug=sXaug, which="LevMar_step", LM_z=zInfo, damping=damping)
        v_h_beta <- list(
          v_h=old_v_h_beta$v_h + LevMarblob$dv_h,
          beta_eta=old_v_h_beta$beta_eta + LevMarblob$dbeta_eta
        )  
      } else if ( which_LevMar_step=="v") { ## v_h estimation given beta 
        LevMarblob <- get_from_MME(sXaug=sXaug, which="LevMar_step_v_h", LM_z=zInfo, damping=damping)
        v_h_beta <- list(
          v_h=old_v_h_beta$v_h + LevMarblob$dv_h,
          beta_eta=old_v_h_beta$beta_eta 
        )  
      }
    } else { ## joint hlik maximization
      LevMarblob <- get_from_MME(sXaug=sXaug, which="LevMar_step", LM_z=zInfo, damping=damping)
      v_h_beta <- list(
        v_h=old_v_h_beta$v_h + LevMarblob$dv_h,
        beta_eta=old_v_h_beta$beta_eta + LevMarblob$dbeta_eta
      )  
    }
    eta <- off + drop(processed$AUGI0_ZX$X.pv %*% v_h_beta$beta_eta) + drop(ZAL %id*% v_h_beta$v_h) ## length nobs 
    
    newmuetablob <- .muetafn(eta=eta,BinomialDen=processed$BinomialDen,processed=processed) 
    neww.resid <- .calc_w_resid(newmuetablob$GLMweights,phi_est)
    
    if (is.null(etaFix$v_h)) { 
      if (GLMMbool) {
        u_h <- v_h_beta$v_h 
        newwranefblob <- wranefblob ## keep input wranefblob since GLMM and lambda_est not changed
      } else {
        u_h_v_h_from_v_h_args <- c(constant_u_h_v_h_args,list(v_h=v_h_beta$v_h))
        u_h <- do.call(".u_h_v_h_from_v_h",u_h_v_h_from_v_h_args)
        if ( ! is.null(attr(u_h,"v_h"))) { ## second test = if u_h_info$upper.v_h or $lower.v_h non NULL
          v_h_beta$v_h <- attr(u_h,"v_h")
        }
        ## update functions u_h,v_h
        newwranefblob <- do.call(".updateW_ranefS",c(updateW_ranefS_constant_arglist,list(u_h=u_h,v_h=v_h_beta$v_h)))
      } 
    } else newwranefblob <- wranefblob
    
    newsXaug <- do.call(processed$AUGI0_ZX$envir$method, # ie, def_AUGI0_ZX_sparsePrecision
                        c(update_sXaug_constant_arglist,
                          list(w.ranef=newwranefblob$w.ranef, 
                               w.resid=neww.resid)))
    APHLs_args$sXaug <- newsXaug
    APHLs_args$dvdu <- newwranefblob$dvdu
    APHLs_args$u_h <- u_h 
    APHLs_args$mu <- newmuetablob$mu
    newAPHLs <- do.call(".calc_APHLs_from_ZX", APHLs_args)
    if (processed$HL[1L]==1L) { 
      if (which_LevMar_step=="v") { 
        newlik <- newAPHLs[["hlik"]]
        oldlik <- oldAPHLs[["hlik"]]
        if (damping==0L) {
          if (trace) print(paste("IRLS step for v_h, hlik=",newlik))
          break
        }
      } else {
        newlik <- newAPHLs[["p_v"]]
        oldlik <- oldAPHLs[["p_v"]]
        if (damping==0L) {
          if (trace) print(paste("IRLS step for (beta,v_h); p_v=",newlik))
          break
        }
      } 
    } else {
      newlik <- newAPHLs[["hlik"]]
      oldlik <- oldAPHLs[["hlik"]]
      if (damping==0L) {
        if (trace) print(paste("IRLS step, hlik=",newlik))
        break
      }
    }

    gainratio <- (newlik!=-Inf) ## -Inf occurred in binary probit with extreme eta... 
    if (gainratio) { 
      if (processed$HL[1L]==1L) { ## ML fit 
        if (which_LevMar_step=="v_b") { 
          dv_h_beta <- unlist(LevMarblob[c("dv_h","dbeta_eta")])
          summand <- dv_h_beta*(gainratio_grad+ LevMarblob$dampDpD * dv_h_beta) 
        } else if (which_LevMar_step=="v") { ## v_h estimation given beta (FIXME can surely be made more exact)
          summand <- LevMarblob$dv_h*(gainratio_grad[seq_n_u_h]+ LevMarblob$dampDpD * LevMarblob$dv_h) 
        }
      } else { ## joint hlik maximization
        dv_h_beta <- unlist(LevMarblob[c("dv_h","dbeta_eta")])
        summand <- dv_h_beta*(gainratio_grad+ LevMarblob$dampDpD * dv_h_beta) 
      }
      ## The two terms of the summand should be positive. In part. dv_h_beta*gainratio_grad should be positive. 
      ## However, numerical error may lead to <0 or even -Inf
      ## Further, if there are both -Inf and +Inf elements the sum is NaN.
      summand[summand<0] <- 0
      denomGainratio <- sum(summand)
      dlogL <- newlik-oldlik
      conv_logL <- abs(dlogL)/(1+abs(newlik))
      gainratio <- 2*dlogL/denomGainratio ## cf computation in MadsenNT04 below 3.14, but generalized for D' != I ; numerator is reversed for maximization
    } else { ## 2017/10/16 a patch to prevent a stop in this case, but covers up dubious computations (FIXME)
      newlik <- -.Machine$double.xmax
      dlogL <- newlik-oldlik
      conv_logL <- abs(dlogL)/(1+abs(newlik))
    }
    if (trace) print(paste("dampingfactor=",dampingfactor,#"innerj=",innerj,
                           "damping=",damping,"gainratio=",gainratio,"oldlik=",oldlik,"newlik=",newlik))
    if (is.nan(gainratio)) {
      # if the initial logL is the solution logL, then damping diverges 
      # it is then possible that some element of dVscaled_beta =0 and some of dampDpD =Inf
      # then the summand has some NaN
      # At the same time not all elements of dVscaled_beta need be 0 (eg small changes in eta for mu=0 or 1 in binomial models)
      # so testing dVscaled_beta is not sufficient to stop the algo
      # (LevenbergM is quite slow in such cases)
      break
    }
    if (gainratio > 0) { 
      ## cf Madsen-Nielsen-Tingleff again, and as in levmar library by Lourakis
      damping <- damping * max(1/3,1-(2*gainratio-1)^3)  
      break 
    } else if (dampingfactor>4 ## implies iter>2
               && conv_logL <1e-8 && abs(prev_conv_logL) <1e-8) { ## FIXME: use less strict criterion ?
      damping <- initdamping ## bc presumably damping has diverged un-usefully
      break ## apparently flat likelihood; this has occurred when fit_as_ZX used a wrong initial (beta,v_h),
      ## but also occurs bc the grad is not that of a single function
    } else { ## UNsuccessful step
      prev_conv_logL <- conv_logL
      damping <- damping*dampingfactor
      dampingfactor <- dampingfactor*2
    }
    if (damping>1e100) stop("reached damping=1e100")
  } 
  return(list(lik=newlik,APHLs=newAPHLs, damping=damping,v_h_beta=v_h_beta, sXaug=newsXaug,
              # fitted=fitted, ## FIXME: removed so that no shortcut a la Bates in calc_APHLs_from_ZX; reimplement the shorcut in that fn?
              eta=eta,muetablob=newmuetablob,wranefblob=newwranefblob,
              u_h=u_h, w.resid=neww.resid))
  
}



.solve_IRLS_as_spprec <- function(ZAL, 
                                   y, n_u_h, 
                       H_global_scale=1, 
                       lambda_est, muetablob=NULL,
                       off, maxit.mean, etaFix,
                       ## for ! LMM
                       eta=NULL, 
                       ## supplement for LevenbergM
                       beta_eta,
                       ## supplement for ! GLMM
                       wranefblob, u_h, v_h, w.resid=NULL, phi_est,
                       for_init_z_args,
                       for_intervals,
                       ##
                       processed, corrPars,
                       trace=FALSE) # corrPars needed together with adjMatrix to define Qmat
{
  pforpv <- ncol(processed$AUGI0_ZX$X.pv)
  nobs <- length(y)
  seq_n_u_h <- seq_len(n_u_h)
  ypos <- n_u_h+seq_len(nobs)
  lcrandfamfam <- attr(processed$rand.families,"lcrandfamfam")
  LMMbool <- processed$LMMbool
  GLMMbool <- processed$GLMMbool
  LevenbergM <- (.determine_LevenbergM(processed$LevenbergM) && is.null(for_intervals))
  is_HL1_1 <- (processed$HL[1L]==1L)
  which_LevMar_step <- "v_b"
  old_relV_beta <- NULL
  not_moving <- FALSE
  damped_WLS_blob <- NULL
  Ftol_LM <- processed$spaMM_tol$Ftol_LM
  if ( LevenbergM) { 
    damping <- 1e-7
    loc_Xtol_rel <- 1e-03 ## maybe good compromise between optim accuracy and time. 
  } else {
    damping <- 0L ## indicator for early exit from .do_damped_WLS without a check for increase 
    loc_Xtol_rel <- processed$spaMM_tol$Xtol_rel/10
  }
  if ( ! LMMbool) {
    constant_zAug_args <- list(n_u_h=n_u_h, nobs=nobs, pforpv=pforpv, y=y, off=off, ZAL=ZAL, processed=processed)
    if ( ! GLMMbool) {
      constant_init_z_args <- c(list(lcrandfamfam=lcrandfamfam, nobs=nobs, lambda_est=lambda_est, ZAL=ZAL),  
                                # fit_as_ZX args specific for ! GLMM:
                                for_init_z_args,
                                #
                                mget(c("cum_n_u_h","rand.families","stop.on.error"),envir=processed))
      constant_u_h_v_h_args <- c(mget(c("cum_n_u_h","rand.families"),envir=processed),
                                 processed$u_h_info, ## elements of u_h_info as elements of constant_u_h_v_h_args  
                                 list(lcrandfamfam=lcrandfamfam))
      updateW_ranefS_constant_arglist <- c(mget(c("cum_n_u_h","rand.families"),envir=processed),list(lambda=lambda_est))
    } 
  } 
  
  ##### initial sXaug
  if (is.null(eta)) { ## NULL input eta allows NULL input muetablob
    eta <- off + drop(processed$AUGI0_ZX$X.pv %*% beta_eta) + drop(ZAL %id*% v_h) 
    muetablob <- .muetafn(eta=eta,BinomialDen=processed$BinomialDen,processed=processed) 
  }
  if ( ! is.null(for_intervals)) {
    v_h_beta <- list(v_h=v_h, beta_eta=for_intervals$beta_eta)
  } else {
    v_h_beta <- list(v_h=v_h,beta_eta=beta_eta)
  }
  ## varies within loop if ! LMM since at least the GLMweights in w.resid change
  if ( is.null(w.resid) ) w.resid <- .calc_w_resid(muetablob$GLMweights,phi_est)
  ## needs adjMatrix and corrPars to define Qmat
  update_sXaug_constant_arglist <- list(AUGI0_ZX=processed$AUGI0_ZX, corrPars=corrPars, 
                                        cum_n_u_h=processed$cum_n_u_h) 
  sXaug <- do.call(processed$AUGI0_ZX$envir$method, # ie, def_AUGI0_ZX_sparsePrecision
                   c(update_sXaug_constant_arglist,
                     list(w.ranef=wranefblob$w.ranef, 
                          w.resid=w.resid)))
  allow_LM_restart <- ( ! LMMbool && ! LevenbergM && is.null(for_intervals) && is.na(processed$LevenbergM["user_LM"]) )
  if (allow_LM_restart) {
    keep_init <- new.env()
    #names_keep <- ls()  
    names_keep <- c("sXaug","wranefblob","muetablob","u_h","w.resid","eta","v_h","beta_eta","v_h_beta")
    for (st in names_keep) keep_init[[st]] <- environment()[[st]]
  }
  LMcond <- - 10. 
  ################ L O O P ##############
  for (innerj in 1:maxit.mean) { ## FIXME gradient => potential ? => single objective function ? 
    if( ! LevenbergM && allow_LM_restart) { ## FIXME the next step improvement would be 
      #  ./. to keep track of lowest lambda that created problem and use LM by default then
      if (innerj>3) {
        LMcond <- LMcond + mean(abs_d_relV_beta/(old_abs_d_relV_beta+1e-8))
        ## cat(mean(abs_d_relV_beta/old_abs_d_relV_beta)," ")
        # cat(LMcond/innerj," ")
        if (LMcond/innerj>0.5) {
          if (trace) cat("!LM")
          for (st in names_keep) assign(st,keep_init[[st]])
          # v_h_beta included in keep_init
          LevenbergM <- TRUE
          damping <- 1e-7
          loc_Xtol_rel <- 1e-03 ## maybe good compromise between optim accuracy and time. 
          damped_WLS_blob <- NULL
          allow_LM_restart <- FALSE 
        }
      }
      if (innerj>2) old_abs_d_relV_beta <- abs_d_relV_beta
    }
    ##### get the lik of the current state
    if ( ! is.null(for_intervals)) {
      loc_logLik_args <- list(sXaug=sXaug, processed=processed, phi_est=phi_est,
                              lambda_est=lambda_est, dvdu=wranefblob$dvdu, u_h=u_h, mu=muetablob$mu)
      oldlik <- unlist(do.call(".calc_APHLs_from_ZX",loc_logLik_args)[for_intervals$likfn]) # unlist keeps name
    } else if (LevenbergM) { ## then logL is necessary to check for increase
      if (is.null(damped_WLS_blob)) {
        oldAPHLs <- .calc_APHLs_from_ZX(sXaug=sXaug, processed=processed, phi_est=phi_est, which=processed$p_v_obj, 
                                        lambda_est=lambda_est, dvdu=wranefblob$dvdu, u_h=u_h, mu=muetablob$mu)
      } else { ## Levenberg and innerj>1
        oldAPHLs <- damped_WLS_blob$APHLs
      }
    } 
    if (LMMbool) {
      zInfo <- list(z2=NULL,z1=y-off,sscaled=0)
      zInfo$z1_eta <- z1_sscaled_eta <- y-eta
    } else {
      if ( ! GLMMbool) {
        # arguments for init_resp_z_corrections_new called in calc_zAug_not_LMM
        init_z_args <- c(constant_init_z_args,
                         list(w.ranef=wranefblob$w.ranef, u_h=u_h, v_h=v_h_beta$v_h, dvdu=wranefblob$dvdu, 
                              sXaug=sXaug, w.resid=w.resid))
      } else init_z_args <- NULL
      calc_zAug_args <- c(constant_zAug_args,
                          list(eta=eta, muetablob=muetablob, dlogWran_dv_h=wranefblob$dlogWran_dv_h, 
                               sXaug=sXaug, ## .get_qr may be called for Pdiag calculation
                               w.ranef=wranefblob$w.ranef, 
                               w.resid=w.resid,
                               init_z_args=init_z_args) )
      zInfo <- do.call(".calc_zAug_not_LMM",calc_zAug_args)
      if (GLMMbool) zInfo$z2 <- NULL
      etamo <- eta - off
      zInfo$z1_eta <- zInfo$z1- etamo 
      z1_sscaled_eta <- zInfo$z1_sscaled - etamo # augz[-seq_n_u_h]- etamo # z_1-sscaled-etamo
    }
    ## keep name 'w'zAug to emphasize the distinct weightings  of zaug and Xaug (should have been so everywhere)
    #####

    ##### improved  Vscaled_beta   
    ## he solver uses a 'beta first" approach to solveing for d_beta and d_v... even for LMM 
    if (GLMMbool) {
      zInfo$dlogfvdv <-  - v_h_beta$v_h * wranefblob$w.ranef
    } else zInfo$dlogfvdv <- (zInfo$z2 - v_h_beta$v_h) * wranefblob$w.ranef
    ## $gainratio_grad is the rhs in the direct solution of the full system (by chol2inv in the dense QR case)
    ## and thus propto the gradient of the objective.
    ## The 'betaFirst' algo in sparsePrecision does not use it, only the gainratio code does  
    ## the gradient for -p_v (or -h), independent of the scaling
    zInfo$m_grad_obj <- c( ## drop() avoids c(Matrix..); attr(sXaug,"w.resid") correct for truncated models too.
      m_grad_v <- drop(.crossprod(ZAL, attr(sXaug,"w.resid") * zInfo$z1_eta) +zInfo$dlogfvdv), # Z'W(z_1-eta) + dlogfvfv
      drop(.crossprod(processed$AUGI0_ZX$X.pv, attr(sXaug,"w.resid") * z1_sscaled_eta)) # X'W(z_1-sscaled-eta) 
    )
    #
    if ( ! is.null(for_intervals)) {
      currentDy <- (for_intervals$fitlik-oldlik)
      if (currentDy < -1e-4) .warn_intervalStep(oldlik,for_intervals)
      intervalBlob <- .intervalStep_spprec(old_v_h_beta=v_h_beta,
                                       sXaug=sXaug,zInfo=zInfo,
                                       for_intervals=for_intervals,
                                       currentlik=oldlik,currentDy=currentDy)
      damped_WLS_blob <- NULL
      v_h_beta <- intervalBlob$v_h_beta
    } else { ## IRLS or not !
      # gradient for scaled system from gradient of objective
      #zInfo$v_h_beta <- v_h_beta
      if (LevenbergM) { ## not IRLS
        # FIXME I made not_moving a sufficent condition fro break below !
        if (not_moving && is_HL1_1) { ## not_moving TRUE may occur when we are out of solution space. Hence test Mg_solve_g
          Mg_solve_g <- sum(zInfo$m_grad_obj*unlist(get_from_MME(sXaug,szAug=zInfo))) ## FIXME presumably not efficient 
          if (Mg_solve_g < Ftol_LM/2) {
            if (trace>1L) {"break bc Mg_solve_g<1e-6"}
            break
          }
        } ## else not_moving was a break condition elsewhere in code
        
        if (trace>1L )  {
          maxs_grad <- c(max(abs(zInfo$m_grad_obj[seq_n_u_h])),max(abs(zInfo$m_grad_obj[-seq_n_u_h])))
          if (trace) cat("iter=",innerj,", max(|grad|): v=",maxs_grad[1L],"beta=",maxs_grad[2L],";")
        }
        constant_APHLs_args <- list(processed=processed, which=processed$p_v_obj, sXaug=sXaug, phi_est=phi_est, lambda_est=lambda_est)
        if (is_HL1_1 && ! is.null(damped_WLS_blob)) {
          if (which_LevMar_step=="v") {
            hlik_stuck <- (damped_WLS_blob$APHLs$hlik < oldAPHLs$hlik + Ftol_LM/10)
            if ( ! hlik_stuck) need_v_step <- (get_from_MME(sXaug=damped_WLS_blob$sXaug, which="Mg_invH_g", B=m_grad_v) > Ftol_LM/2) 
            if ( hlik_stuck || ! need_v_step) { ## LevMar apparently maximized h wrt v after several iterations
              ## if hlik has not recently moved or has moved but reached a point where the h gradient is negligible 
              if (trace>2L) print("switch from v to v_b")
              old_relV_beta <- relV_beta ## serves to assess convergence !!!!!!!!!!!!!!!!!!!
              which_LevMar_step <- "v_b" 
            } else {
              if (trace>2L) print("still v")
              ## v_h estimation not yet converged, continue with it
            }
          } else { ## performed one improvement of p_v by new v_b, 
            # indirect way of checking Mg_solve_g:
            # FIXME I made not_moving a sufficent condition fro break below !
            if (not_moving) { ## if we reach this point, Mg_solve_g (tested above) was too large, we must be out of solution space
              # need_v_step <- TRUE ## implicit meaning
            } else {
              p_v_stuck <- (damped_WLS_blob$APHLs$p_v < oldAPHLs$p_v + Ftol_LM/10) ## test whether LevMar apparently solved (v,beta) equations after several iterations
              if ( ! p_v_stuck) need_v_step <- (get_from_MME(sXaug=damped_WLS_blob$sXaug, which="Mg_invH_g", B=m_grad_v) > Ftol_LM/2) 
              ## we have identified two gradient cases where we must check v: Mg_solve_g>0 or (if estimates have just moved) Mg_invH_g>0 
            }
            if ( not_moving || p_v_stuck || need_v_step) { ## logically we may not need p_v_stuck, but this condition is faster to evaluate
              # p_v_stuck is analogous to (not_moving BUT large Mg_solve_g), checking that lik and estimates do not change 
              if (trace>2L) print("switch from v_b to v")
              which_LevMar_step <- "v"
            } else {
              if (trace>2L) print("still v_b")
            }
          }
        }         
        damped_WLS_blob <- .do_damped_WLS_spprec(sXaug=sXaug, 
                                                 zInfo=zInfo, 
                                                 old_v_h_beta=v_h_beta, old_eta=eta,
                                                 oldAPHLs=oldAPHLs,
                                                 APHLs_args = constant_APHLs_args,
                                                 damping=damping,
                                                 ypos=ypos,off=off,
                                                 GLMMbool=GLMMbool,etaFix=etaFix,
                                                 constant_u_h_v_h_args=constant_u_h_v_h_args,
                                                 updateW_ranefS_constant_arglist=updateW_ranefS_constant_arglist,
                                                 wranefblob=wranefblob,seq_n_u_h=seq_n_u_h,
                                                 update_sXaug_constant_arglist=update_sXaug_constant_arglist, 
                                                 processed=processed, 
                                                 phi_est=phi_est, n_u_h=n_u_h, ZAL=ZAL,
                                                 which_LevMar_step = which_LevMar_step
        )
        ## LevM PQL
        if (! is_HL1_1) {
          if (damped_WLS_blob$lik < oldAPHLs$hlik) { ## if LevM step failed to find a damping that increases the lik
            ## This occurs inconspiscuously in the PQL_prefit providing a bad starting point for ML fit
            damped_WLS_blob <- NULL
            dv_h_beta <- get_from_MME(sXaug,szAug=zInfo) ################### FIT
            v_h_beta <- list(v_h=v_h_beta$v_h+dv_h_beta$dv_h,
                             beta_eta=v_h_beta$beta_eta+dv_h_beta$dbeta_eta)
            LevenbergM <- FALSE ## D O N O T set it to TRUE again !
          } 
        }
      } else { ## IRLS: always accept new v_h_beta
        damped_WLS_blob <- NULL
        dv_h_beta <- get_from_MME(sXaug,szAug=zInfo) ################### FIT
        v_h_beta <- list(v_h=v_h_beta$v_h+dv_h_beta$dv_h,
                      beta_eta=v_h_beta$beta_eta+dv_h_beta$dbeta_eta)
      }
    }
    ######
    
    ##### Everything that is needed for 
    #  (1) assessment of convergence: c(v_h*sqrt(wranefblob$w.ranef),beta_eta)
    #  (2) all return elements are updated as function of the latest Vscaled_beta.
    #      In particular We need muetablob and (if ! LMM) sXaug, hence a lot of stuff.
    #  Hence, the following code is useful whether a break occurs or not. 
    if ( ! is.null(damped_WLS_blob) ) {
      v_h_beta <- damped_WLS_blob$v_h_beta 
      eta <- damped_WLS_blob$eta
      wranefblob <- damped_WLS_blob$wranefblob
      u_h <- damped_WLS_blob$u_h
      muetablob <- damped_WLS_blob$muetablob
      w.resid <- damped_WLS_blob$w.resid
      sXaug <- damped_WLS_blob$sXaug
    } else {
      # drop, not as.vector(): names are then those of (final) eta and mu -> used by predict() when no new data
      eta <- off + drop(processed$AUGI0_ZX$X.pv %*% v_h_beta$beta_eta) + drop(ZAL %id*% v_h_beta$v_h) ## length nobs 
      if ( is.null(etaFix$v_h)) {
        if (GLMMbool) {
          u_h <- v_h_beta$v_h ## keep input wranefblob since lambda_est not changed
        } else {
          u_h_v_h_from_v_h_args <- c(constant_u_h_v_h_args,list(v_h=v_h_beta$v_h))
          u_h <- do.call(".u_h_v_h_from_v_h",u_h_v_h_from_v_h_args)
          if ( ! is.null(attr(u_h,"v_h"))) { ## second test = if u_h_info$upper.v_h or $lower.v_h non NULL
            v_h_beta$v_h <- attr(u_h,"v_h")
          }
          ## update functions u_h,v_h
          wranefblob <- do.call(".updateW_ranefS",c(updateW_ranefS_constant_arglist,list(u_h=u_h,v_h=v_h_beta$v_h)))
        }
      }
      muetablob <- .muetafn(eta=eta,BinomialDen=processed$BinomialDen,processed=processed) 
      if ( ! LMMbool ) {
        w.resid <- .calc_w_resid(muetablob$GLMweights,phi_est)
        sXaug <- do.call(processed$AUGI0_ZX$envir$method, # ie, def_AUGI0_ZX_sparsePrecision
                         c(update_sXaug_constant_arglist, 
                           list(w.ranef=wranefblob$w.ranef, 
                                w.resid=w.resid)))
      }
    }
    #####
    
    ##### assessment of convergence
    if (innerj<maxit.mean) {
      relV_beta <- c(v_h_beta$v_h*sqrt(wranefblob$w.ranef),v_h_beta$beta_eta)  ## convergence on v_h relative to sqrt(lambda), more exactly for Gaussian
      abs_d_relV_beta <- abs(relV_beta - old_relV_beta)
      not_moving <- ( ( ! is.null(old_relV_beta)) && mean(abs_d_relV_beta) < loc_Xtol_rel )
      if (is.na(not_moving)) {
        if (anyNA(relV_beta)) {
          if ( ! is.null(damped_WLS_blob)) {
            message(paste("innerj=",innerj,"damping=",damping,"lik=",damped_WLS_blob$lik))
            stop("Numerical problem despite Levenberg algo being used: complain.")
          } else stop("Numerical problem: try control.HLfit=list(LevenbergM=TRUE)")
        } else stop("Error in evaluating break condition")
      } 
      if (not_moving) break ## sufficient condition here
      if ( ! (is_HL1_1 && LevenbergM)) { ## possible reversal of condition from F to T in  LevM PQL !!!!
        old_relV_beta <- relV_beta
      } ## ELSE old_relV_beta controlled in block for which_LevMar_step !!
    } else break
  } ################ E N D LOOP ##############
  
  names(v_h_beta$beta_eta) <- colnames(processed$AUGI0_ZX$X.pv)
  RESU <- list(sXaug=sXaug, 
               ## used by calc_APHLs_from_ZX
               #fitted=fitted, ## FIXME: removed so that no shortcut a la Bates in calc_APHLs_from_ZX; reimplement the shorcut in that fn?  
               weight_X=NA, nobs=nobs, pforpv=pforpv, seq_n_u_h=seq_n_u_h, u_h=u_h, 
               muetablob=muetablob, 
               lambda_est=lambda_est,
               phi_est=phi_est,
               ## used by other code
               beta_eta=v_h_beta$beta_eta, w.resid=w.resid, wranefblob=wranefblob, 
               v_h=v_h_beta$v_h, eta=eta, innerj=innerj)
  return(RESU)
} 

.intervalStep_spprec <- function(old_v_h_beta,sXaug,zInfo,currentlik,for_intervals,currentDy) {
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
    ## (fixme) patch:
    if (currentDy<targetDy) { ## we are close to the ML: we extrapolate a bit more confidently
      min_abs_Dx <- for_intervals$asympto_abs_Dparm/1000
    } else min_abs_Dx <- 1e-6 ## far from ML: more cautious move our of Dx=0
    Dx <- sign(currentDx)*max(abs(Dx),min_abs_Dx)
    v_h_beta_vec[parmcol_ZX] <- for_intervals$MLparm + Dx 
  }
  ## szAug changes:
  parmcol_X <- for_intervals$parmcol_X
  zInfo$z1 <- zInfo$z1 - sXaug$AUGI0_ZX$X.pv[,parmcol_X]*v_h_beta_vec[parmcol_ZX] 
  zInfo$m_grad_obj <- zInfo$m_grad_obj[-parmcol_ZX]
  ## sXaug changes: 
  int_AUGI0_ZX <- as.list(sXaug$AUGI0_ZX) ## as.list() local copy avoids global modifs of original sXaug$AUGI0_ZX envir
  int_AUGI0_ZX$X.pv <- int_AUGI0_ZX$X.pv[,-(parmcol_X),drop=FALSE]
  int_AUGI0_ZX$ZeroBlock <- int_AUGI0_ZX$ZeroBlock[,-(parmcol_X),drop=FALSE]
  # Assuming sXaug not being an envir, sXaug$AUGI0_ZX <- int_AUGI0_ZX would be correct, creating a local copy of sXaug, 
  #  but less explicit than the following: 
  int_sXaug <- sXaug ## Assuming sXaug not being an envir, then the following does not affect sXaug.
  int_sXaug$AUGI0_ZX <- int_AUGI0_ZX ## Replaces an envir by a list i nthe local copy; 
  v_h_beta_vec[-(parmcol_ZX)] <- old_v_h_beta_vec[-(parmcol_ZX)]+ unlist(get_from_MME(int_sXaug,szAug=zInfo)) 
  return(list(v_h_beta=relist(v_h_beta_vec,old_v_h_beta))) 
}

