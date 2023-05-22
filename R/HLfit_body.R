HLfit_body <- function(processed, 
         control.HLfit=list(), ## used both by preprocess and HLfit_body
         init.HLfit = list(), ## not from processed: this is affected by HLCor_body -> .canonizeRanPars(ranPars) post .preprocess()ing
         #                       Thus if in a HLCor call we can expect a $corrPars in sp3 code.
         fixed=list(), ##  possibly trLambda but not necessarily.
         etaFix=list() ## beta, v_h (or even u_h)
) {
  processed$envir$ranFix <- fixed # for diagnostics reported by div_info() (seek '$ranFixes') [_F I X M E_ rethink] 
  ranFix <- .post_process_respfamilies(processed$family, ranFix=fixed, families=processed$families) ## assign 'extra' COMPoisson or negbin pars and cleans ranFix of them
  # next line to be called before we extract anything (lambda, ranCoefs... ) from ranFix:
  ranFix <- .canonizeRanPars(ranPars=ranFix,corr_info=NULL, checkComplete = FALSE, rC_transf=.spaMM.data$options$rC_transf)## including full-size lambda
  #data <- processed$data
  verbose <- processed$verbose
  
  predictor <- attr(processed$predictor,"no_offset") 
  # prior.weights <- processed$prior.weights
  
  warningEnv <- new.env(parent=emptyenv())
  ## when adding verbose elements, remind that these might be lost through corrHLfit -> HLCor cf dotlist$verbose <- verbose[intersect(...]
  ##
  y <- processed$y
  HL <- processed$HL
  spaMM_tol <- processed$spaMM_tol
  max.iter <- processed$max.iter
  BinomialDen <- processed$BinomialDen
  #X.pv <- processed$AUGI0_ZX$X.pv
  ### a bit of post processing
  pforpv <- ncol(processed$AUGI0_ZX$X.pv)
  models <- processed$models
  LMMbool <- attr(models,"LMMbool")
  #### Note that HLCor modifies the L matrix (=> ZAL cannot be preprocessed by corrHLfit and must be recomputed each time 
  ZAlist <- processed$ZAlist ## : ZAlist is a list of design matrices 
  nrand <- length(ZAlist)
  cum_n_u_h <- processed$cum_n_u_h
  n_u_h <- cum_n_u_h[nrand+1L] 
  vec_n_u_h <- diff(cum_n_u_h)
  if (nrand) { 
    ranCoefs_blob <- .process_ranCoefs(processed, ranCoefs=.getPar(ranFix,"ranCoefs"), ## may be NULL, 
                                       use_tri_CORREL=TRUE) ## UPDATES preexisting object # no augZXy pb here
    LMatrices <- .process_LMatrices(processed, ranCoefs_blob)
    which_inner_ranCoefs <- attr(LMatrices,"which_inner_ranCoefs") # needed & non-NULL only for spprec
    ZAL <- .process_ZAL(processed, LMatrices, ZAlist, HL)
  } else ZAL <- NULL 
  ## 
  ranFix$lambda <- .reformat_lambda(ranFix$lambda, nrand, namesTerms=attr(ZAlist,"namesTerms"), full_lambda=TRUE) # necessary to standardize names before next line
  if (any(ranFix$lambda==0,na.rm=TRUE)) stop("lambda cannot be fixed to 0.")
  lam_fix_or_outer_or_NA <- processed$reserve$repNAnrand
  lam_fix_or_outer_or_NA[names(ranFix$lambda)] <- ranFix$lambda # .getPar(ranFix,"lambda") ## should already have length 'nrand' or else be NULL
  ###################
  if (! is.null(corr_est <- .get_cP_stuff(init.HLfit,"rho"))) corr_est <- list(rho=corr_est) 
  if (need_simple_lambda <- need_ranefPars_estim <- (nrand>0L)) {
    need_simple_lambda <- any(is.na(lam_fix_or_outer_or_NA) & ! ranCoefs_blob$is_set)
    need_ranefPars_estim <-  (need_simple_lambda || ! is.null(corr_est))
  } 
  # ranFix overrides $phi.Fix so that the HLCorcall can be used in post-fit code to compute numerical info matrix:
  # but... the constr_phi and constr_fit attributes of $phi.Fix are lost, while constr_fit (at least) is needed in post-fit
  # So they must be put back by .add_phi_returns() -> .get_phi_object(() in a full fit object.
  if (is.list(phi.Fix <- processed$phi.Fix)) {
    phi.Fix <- .modify_list(phi.Fix,.getPar(ranFix,"phi"))
  } else if (is.null(phi.Fix <- .getPar(ranFix,"phi"))) phi.Fix <- processed$phi.Fix
  ## => initial value is preprocessed value. If the latter is NULL, this remains NULL, 
  # except when  RHS was set in final call of outer estimation 
  # (it would be misleading to compute leverages in such a final call)
  # ..modify_list() is necessary for fitmv: phi.Fix <- .getPar(ranFix,"phi") may replace a partially NULL by a full NULL.
  # This however means that (internally-set) ranFix phi component should be a *named list*.
  
  nothing_to_fit <-  ((! need_ranefPars_estim) && pforpv==0L && (! .anyNULL(phi.Fix)) 
                      && (nrand && (! is.null(etaFix$v_h))) )
  if ( nothing_to_fit ) { 
    whichadj <- which(attr(ZAlist,"exp_ranef_types")=="adjacency") ## bug presumably corrected here 30/12/2017
    fixed_adjacency_info <- .get_fixed_adjacency_info(whichadj, LMatrices, cum_n_u_h, corr_est, ranFix, init.HLfit)
    # only APHLs:
    return(.nothing_to_fit(phi.Fix, off, models, etaFix, processed$rand.families, cum_n_u_h, 
                           lam_fix_or_outer_or_NA, vec_n_u_h, n_u_h, fixed_adjacency_info, ZAL, BinomialDen, processed)) 
    # => Possible error with .do_TRACE bc the exit tracing code does not find the 'res' variable, not locally defined in the case. I could add res <- ... here.
  }   ### RETURN !! ## not of class HLfit, and p_bv is not returned.
  
  ##########################################################################################
  ##########################################################################################
  ##########################################################################################
  std_dev_res_needed_4_inner_estim <- .anyNULL(phi.Fix) || need_simple_lambda 
  ### Initial values  for lambda, phi and beta from lam_fix_or_outer_or_NA, phi.Fix, or init.HLfit ##### 
  ## Initial estimate for phi 
  vec_nobs <- processed$vec_nobs
  phi_est <- phi.Fix
  phi_est <- .denullify(phi_est, modifier=processed$port_env$port_fit_values$phi_est, vec_nobs=vec_nobs)
  phi_est <- .denullify(phi_est, modifier=init.HLfit$phi, vec_nobs=vec_nobs)
  if (nrand) { 
    if (identical(processed$return_only,"p_vAPHLs")) {
      whichAPHLs <- "p_v"
    } else if (identical(processed$return_only,"p_bvAPHLs")) {
      whichAPHLs <- "p_bv" ## return value may still include p_v if it is used to compute p_bv
    } else whichAPHLs <- c("p_v","p_bv")
    ## Initial estimate for lambda in 'compact" form:
    init.lambda <- .calc_initial_init_lambda(lam_fix_or_outer_or_NA, nrand, processed, ranCoefs_blob, init.HLfit, ranFix)
  } else {
    # u_h <- v_h <- numeric(0)
    init.lambda <- NULL
  }
  ###
  ## Initial estimate for beta  (etaFix does NOT act directly in .wrap_IRLS -> .solve_IRLS...)
  ###
  if ( ! is.null(processed$X_off_fn)) { # (___F I X M E___?) currently X_off_fn does not allow partial beta's (with potential mess with initial beta_eta )
    beta_eta <- numeric(0)
    processed$off <- off <- processed$X_off_fn(etaFix$beta) # .solve_IRLS_as_ZX() uses processed$off
    # AUGI0_ZX$X.pv must correspondly have been reduced by .preprocess
  } else {
    off <- processed$off
    beta_eta <- .get_init_beta(processed, pforpv, init.HLfit) # (note that this correctly avoids is.null(beta_eta) ***when*** pforpv=0) 
                                                              # __F I X M E__ what do we exactly need for LMMs (?) 
  }
  ######### missing Initial estimates for mu, phi, lambda by GLM ####################
  if ( is.null(beta_eta) ||  # occurs when pforpv>0 and .get_init_beta() did not find anything
       .anyNULL(phi_est) || anyNA(init.lambda) ) { 
    inits_by_xLM <- .get_inits_by_xLM(processed, 
                                      reset=quote(family$family %in% c("COMPoisson","negbin1","negbin2", "beta_prec")) ) # quoted to be applied to each family in mv case
    ## : uses processed$y, $BinomialDen, [["control.glm"]]
  }
  if (is.null(beta_eta) ) beta_eta <- inits_by_xLM$beta_eta # from .lm.fit or lm.fit using scaled X.pv, hence result is scaled value.
  #
  if (HL[1]=="SEM") {
    
    loopout_blob <-  .wrap_wrap_SEM(processed, ZAL, beta_eta, off, corr_est, init.lambda, lam_fix_or_outer_or_NA, 
                                                LMatrices, verbose, BinomialDen, phi_est)
    SEMblob <- loopout_blob$SEMblob
    APHLs <- loopout_blob$APHLs
  } else {
    
    intervalInfo <- processed$intervalInfo
    if (!is.null(intervalInfo)) {
      parmcol <- attr(intervalInfo$parm,"col")
      beta_eta[parmcol] <- intervalInfo$init ## already appropriately scaled if X.pv has been scaled
    }  
    ###
    ## predictor from initial values. When there an etaFix, it is here in the 'off'set, and the dims of X.pv and beta_eta are here correspondingly reduced.
    if (nrand) { # (models[["eta"]]=="etaHGLM") {
      ## Initial estimate for u_h, v_h 
      v_h <- intervalInfo$init_v_h
      if (is.null(v_h)) v_h <- .initialize_v_h(processed, etaFix=etaFix, init.HLfit=init.HLfit) ## checks init.HLfit$v_h
      eta <- off + drop(processed$AUGI0_ZX$X.pv %*% beta_eta) + drop(ZAL %id*% v_h)
    } else  eta <- off +drop(processed$AUGI0_ZX$X.pv %*% beta_eta) ## no iteration hence no updating  ## FREQS
    ## conversion to mean of response variable (COUNTS for binomial)
    
    if ( ! is.null(intervalInfo)) {
      if (nrand) { # (models[["eta"]]=="etaHGLM") {
        intervalInfo$ranFix <- ranFix
      } else if (models[[1L]] %in% c("etaGLM")) {
        intervalInfo$parmcol_X <- parmcol 
      }
    }
    maxit.mean <- .maxit.mean(nrand, pforpv, etaFix, LMMbool, intervalInfo, 
                              # for inner ranCoef estim [eg, HLfit3 example with family=Gammalog)] no leverages needed => std_dev_res_needed_4_inner_estim is FALSE 
                              # yet if there are inner_ranCoefs and only inner estim, not outer, is used, iter_mean_dispVar seems appropriate. 
                              # => second condition on which_inner_ranCoefs: [not a strict check of no outer estime but probably equivalent in practice]
                              inner_est_disp_pars=std_dev_res_needed_4_inner_estim || length(which_inner_ranCoefs), 
                              processed, models, phi.Fix)
    
    ## Create the environment that will contain the results of the loop, storing only variables that are (possibly) modified by it.
    loopout_blob <- new.env(parent=emptyenv())
    # Variables used before and after the loop:
    if (nrand) {
      loopout_blob$LMatrices <- LMatrices
      loopout_blob$v_h <- v_h
    }
    loopout_blob$corr_est <- corr_est
    loopout_blob$beta_eta <- beta_eta
    # Add variables not used before the loop:
    .add_loopout_vars(loopout_blob, phi_est, inits_by_xLM, vec_nobs, eta, BinomialDen, processed, nrand, models, 
                      init.lambda, ZAL, cum_n_u_h, vec_n_u_h, n_u_h, ranCoefs_blob, etaFix, init.HLfit)

    ########################################
    ######### Main loop ####################
    ########################################
    #  The loop: modifies the loopout_blob env:
    .loop_while_TRUE(processed=processed, 
                     loopout_blob=loopout_blob, # changes in loopout_blob accessible afterwards 
                     warningEnv=warningEnv,
                     #
                     # variables that are not needed after the loop possibly modified them:
                     nrand=nrand,
                     n_u_h=n_u_h,
                     intervalInfo=intervalInfo,
                     pforpv=pforpv,
                     maxit.mean=maxit.mean,
                     etaFix=etaFix,
                     ranFix=ranFix,
                     phi.Fix=phi.Fix,
                     init.HLfit=init.HLfit,
                     control.HLfit=control.HLfit,
                     verbose=verbose,
                     spaMM_tol=spaMM_tol,
                     need_ranefPars_estim=need_ranefPars_estim,
                     lam_fix_or_outer_or_NA=lam_fix_or_outer_or_NA,
                     need_simple_lambda=need_simple_lambda,
                     LMMbool=LMMbool,
                     max.iter=max.iter,
                     std_dev_res_needed_4_inner_estim=std_dev_res_needed_4_inner_estim,
                     off=off,
                     which_inner_ranCoefs=which_inner_ranCoefs,
                     whichAPHLs=whichAPHLs,
                     ranCoefs_blob=ranCoefs_blob,
                     ZAL=ZAL # possibly modified in the loop but to used afterwards
    ) 

    # trying2avoidlocalvars <- c("APHLs", "conv.phi", "conv.lambda", "conv.corr","conv_logL", "iter", "w.resid",
    #                            "beta_eta", # the loopout copy is scaled, in contrast to the local copy created below 
    #                            "prior.weights", "corr_est", "calcRanefPars_blob", "PHIblob", "u_h", "v_h",
    #                            "auglinmodblob", "leverages", "lambda_est", "phi_est", "LMatrices", "wranefblob",
    #                            #
    #                            "muetablob" # the one copied...
    #                            )
    # if (length(setdiff(ls(loopout_blob), trying2avoidlocalvars))) stop("poorly organized code...")
    ###  for (st in setdiff(ls(loopout_blob), trying2avoidlocalvars)) assign(st, loopout_blob[[st]]) 

  }
  ########################################
  ######### END main loop ################
  ########################################
  if (verbose["trace"]) {
    if ((iter <- loopout_blob$iter)==max.iter) {
      message(paste("(beta,v)/lambda/phi iterations failed to converge in",max.iter,"iterations"))
    } else {
      message(paste("(beta,v)/lambda/phi iterations in HLfit() converged in",iter,"iterations"))
    }
  }
  #
  ######################### potential R E T U R N here: cases without p_bv
  if ( identical(processed$return_only,"p_vAPHLs")) {
    # Following comment no longer clear, but this may have referred to the non-existence of 'processed'  in optimthrousmooth code
    # a bit of ugly coding, but optimthroughsmooth calls HLCor, not HLCor.obj, thus it cannot directly control return_only. So either leave as is, or move the test to HLCor, or modify optimthroughsmooth to call HLCor.obj  
    if (HL[1]=="SEM") { # lambda used for smoothing.
      res <- list(APHLs=APHLs,lambda=SEMblob$lambda) 
    } else {
      APHLs <- loopout_blob$APHLs
      if ( ! is.null(oldp_v <- processed$port_env$objective)) {
        .update_port_fit_values(old_obj=oldp_v,new_obj=APHLs$p_v, loopout_blob=loopout_blob,
                                models=models, processed=processed, control.HLfit=control.HLfit)
      } else if ( ! identical(control.HLfit$write_port_env,FALSE)) {
        processed$port_env$objective <- APHLs$p_v # when return only p_v => not for residModel
      }
    }
    res <- list(APHLs=APHLs) 
    return(res)   ########################   R E T U R N
  } 
  ## ELSE continue: make sure p_bv is included
  if (HL[1] != "SEM") {
    if (nrand) { # (models[["eta"]]=="etaHGLM") {
      APHLs <- loopout_blob$APHLs
      ## cf notes 19/08/2016 pour calcul APHLs et IC's for phiHGLM 
    } else  APHLs <- .calc_APHLs_XLM(processed, w.resid=loopout_blob$w.resid, clik=loopout_blob$APHLs$clik) ## GLM|LLM. 'w.resid' used only for 'REML' which is a new concept for a fixed-effect LLM 
  }
  ######################### potential R E T U R N here: with p_bv
  if ( identical(processed$return_only,"p_bvAPHLs")) {
    if ( ! is.null(oldp_bv <- processed$port_env$objective)) {
      .update_port_fit_values(old_obj=oldp_bv,new_obj=APHLs$p_bv, loopout_blob=loopout_blob,
                              models=models, processed=processed, control.HLfit=control.HLfit)
    } else if ( ! identical(control.HLfit$write_port_env,FALSE)) {
      processed$port_env$objective <- APHLs$p_bv # when return only p_bv => not for residModel
    }
    res <- list(APHLs=APHLs)
    return(res)    ########################   R E T U R N
  } else if (identical(processed$return_only,"confint_bound")) {
    res <- list(APHLs=APHLs)
    res <- .add_unscaled_X.pv_fixef(res=res, processed=processed, beta_eta=loopout_blob$beta_eta, etaFix=etaFix)
    return(res)    ########################   R E T U R N fixef + APHLs
  }
  
  ######################
  ######################
  ######################
  ######################################
  ## BUILD full RETURN VALUE
  ######################################
  ## If an outer optimizer has been called,
  #  "fix" and "outer" parameters are given these types by .get_refit_args() after the optimization call, 
  # then HLfit called again and we reach this point.
  #  This means ranFix gets its type from there *if* properly retained by .canonizeRanPars() 
  #  Then we add inner-optimized parameters, with "var" type added by .get_CorrEst_and_RanFix()
  if ( ! is.null(corr_est) && ! is.null(init.HLfit$corrPars)) corr_est <- list(corrPars=relist(loopout_blob$corr_est$rho,init.HLfit$corrPars)) ## not yet spaMM 3.0
  # Canonical, and inherits all info about outer-optimized corrPars through HLfit's ranFix argument:
  
  CorrEst_and_RanFix <- .get_CorrEst_and_RanFix(ranFix, corr_est) # corr_est parameters are inner-estimated and of type "var"
  how <- list(spaMM.version=packageVersion("spaMM"),
                  MME_method=.get_MME_method(loopout_blob$auglinmodblob, HL=HL),
                  switches=c(augZXy_cond=processed$augZXy_cond, ADFun=processed$ADFun,
                             use_spprec_QR=.spaMM.data$options$use_spprec_QR),
                  obsInfo=processed$how$obsInfo )
  
  res <- list(
    APHLs=APHLs,   ## LIKELIHOODS
    fv=.mu_U2fv(BinomialDen=BinomialDen, muetablob=loopout_blob$muetablob, processed=processed),  ## FITTED VALUES
    muetablob=loopout_blob$muetablob, # directly as a lot of elements may be needed:
    # $dmudeta, eg for dvdlogphiMat in.get_logdispObject() 
    # $p0=muetablob$p0, for simulate 
    # Md3logcLdeta3, Md2logcLdeta2 for hatvalues() -> .hatvals2std_lev() -> . -> ..calc_dlW_deta()
    # $mv for and hatvalues.HLfit(). (and possible for mv simulate)
    eta=.format_eta(eta=loopout_blob$muetablob$sane_eta, data=processed$data), ## convenient for defining starting values... and also sometimes used by predict() # now kept in muetablob
    #
    data=processed$data,
    y=processed$y, 
    ###################
    ## MODEL info      
    ###################
    family=processed$family,
    families=processed$families,
    rand.families=processed$rand.families, 
    models=processed$models, 
    main_terms_info=processed$main_terms_info, ## used by predict
    predictor=processed$predictor, ##  all post fitting functions expect PROCESSED predictor
    vec_nobs=processed$vec_nobs, ## non-null for fitmv
    REMLformula=processed$REMLformula,  # only for .REMLmess()... but it's still a simple way to pass the info. Perhaps put it elsewhere in res?
    distinctX.Re=processed$X.Re, ## NULL if not distinct from X.pv
    #
    QRmethod=processed$QRmethod,
    #
    HL=HL, ## info on fitting objective  # identical to processed$HL
    ranFix=ranFix,
    CorrEst_and_RanFix=CorrEst_and_RanFix,
    how=how)
  
  res <- .add_unscaled_X.pv_fixef(res=res, processed=processed, beta_eta=loopout_blob$beta_eta, etaFix=etaFix)  ## FIXEF, UNSCALED X
  
  if (is.null(processed$family) || processed$family$family %in% c("binomial","betabin")) { # null for mv case
    res$BinomialDen <- BinomialDen # we could put it in all cases...
  }
  
  prior.weights <- processed$prior.weights # with attrs, and possibly a call
  prior.weights[] <- eval(loopout_blob$prior.weights, envir = loopout_blob) ## see Gamma()$simulate # eval the call but keep the attributes.
  res$prior.weights <- prior.weights
  
  .canonize_disp_envs(fitobject=res) # .unscaling resid.models' X
  if ( ! is.null(CorrEst_and_RanFix$corrPars)) {
    res$corrPars <- structure(CorrEst_and_RanFix$corrPars, # ## subset of the above: F I X M E (?) redundancy but convenient when examining fits
                              type=attr(CorrEst_and_RanFix,"type")$corrPars,
                              message='Use get_ranPars(.,which="corrPars") to extract "corrPars" cleanly from fit object.')
  }  
  #
  ##### LAMBDA and other RANEF PARS
  if (need_ranefPars_estim) { # (FALSE if only outer_ranCoefs)
    process_resglm_blob <- .bloc_lambda(processed=processed, lam_fix_or_outer_or_NA=lam_fix_or_outer_or_NA, 
                                        SEMblob=SEMblob, # only for SEM
                                        loopout_blob=loopout_blob, HL=HL)
  } else {
    process_resglm_blob <- list(lambda_pred_list=as.list(rep(NA,nrand)))
  }
  res$dfs <- .calc_dfs(need_ranefPars_estim, process_resglm_blob, init.lambda, 
                   ranCoefs_blob, loopout_blob$LMatrices, processed, pforpv, CorrEst_and_RanFix=CorrEst_and_RanFix)
  
  res$spaMM.version <- structure(res$how$spaMM.version, ## this is NOT a string and comparison with a string is suitably def'ed (as detailed in ?package_version)
                                 message="Please use how(<fit object>)[['spaMM.version']] to extract this information cleanly.")                 
  ###################
  ## LEVERAGES and REML (ie either phi OR lambda was estimated)
  ###################
  if (HL[1L]=="SEM") {
    res$SEM_info <- SEMblob$SEM_info ## info
  } else { ## both lev_phi and deviance_residual missing otherwise
    if (std_dev_res_needed_4_inner_estim) { ## ll model leverages are computed and it makes sense to consider the residuals
      # possible outputs from this block: res$lev_phi; res$lev_lambda; res$std_dev_res; dev_res_blob (also needed out of this block); and info in warningEnv
      res$lev_phi <- loopout_blob$leverages$resid
      dev_res_blob <- .std_dev_resids(res, phi_est=loopout_blob$phi_est, lev_phi=loopout_blob$leverages$resid)  # also needed out of this block
      mu <- res$muetablob$mu
      if (inherits(mu,"Matrix")) {
        warning("inefficiency detected. Please contact the package maintainer.", immediate. = TRUE) # it's inefficient if true in the loop...
        mu <- drop(mu) ## Old comment: "pb calcul deviance_residual" which is why I moved the test from the main loop to here.
      }
      res$std_dev_res <- sign(y-mu) * dev_res_blob$std_dev_res
      
      if (need_simple_lambda) res$lev_lambda <- loopout_blob$leverages$ranef # __F I X M E__ remove the local condition ?
      
      # res$diagnostics$m_grad_obj <- auglinmodblob$m_grad_obj # typically NULL for LMM
      if (nrand && is.null(warningEnv$leveLam1)) { # .calcRanefPars was not called
        # Cf singfitF test in test-rank: outer estim without refit => reaches this point
        #But we want to diagnose possible non-identifiability of lambda => add info to warningEnv:
        .diagnose_lev_lambda(loopout_blob$leverages, warningEnv, nrand, cum_n_u_h=processed$cum_n_u_h)
      }
    }
  }  
  ###################
  ## ALL other LAMBDA returns
  ###################
  ## $w.ranef and $w.resid not doc'ed, as there is no mention of the augmented model in the doc.
  # if (is.list(w.resid)) { ## truncated 'family', or 'families' with some truncated one(s), but not all mv cases
  #   res$w.resid <- w.resid$w_resid ## useful for .get_info_crits() and get_LSmatrix
  # } else res$w.resid <- w.resid ## useful for .get_info_crits() and get_LSmatrix()
  # # res$w.resid is always a vector
  # if (is.null(attr(res$w.resid,"unique"))) attr(res$w.resid,"unique") <- length(unique(res$w.resid))==1L # is.null() => for mv    
  # #
  if (nrand) { # (models[["eta"]]=="etaHGLM") {
    res <- .add_ranef_returns(res, processed=processed, process_resglm_blob=process_resglm_blob, init.lambda=init.lambda,
                              loopout_blob=loopout_blob,
                              ranCoefs_blob=ranCoefs_blob, moreargs=attr(fixed,"moreargs"))
  } ## else various res$ elements are NULL
  ###################
  ## ALL other PHI returns
  ###################
  res <- .add_phi_returns(res=res, processed=processed, loopout_blob=loopout_blob, phi.Fix=phi.Fix, dev_res=dev_res_blob$dev_res)
  ################### the magic environment
  res$envir <- .add_fitobject_envir(nrand=nrand, HL=HL, loopout_blob=loopout_blob, processed=processed)
  ###################
  ## WARNINGS
  ###################
  .hack_options_error(message=NULL)
  ## translation of warnings in user-more friendly form 
  res$warnings <- .post_process_warningEnv(warningEnv=warningEnv, processed=processed, maxit.mean=maxit.mean,  pforpv=pforpv, 
                                           nonSPD=.BLOB(res$envir$sXaug)$nonSPD, 
                                           loopout_blob=loopout_blob # for conv_logL, iter, conv.lambda, conv.phi, conv.corr, ...innerj
  )
  .verbose_warnings(verbose, res$warnings) # immediate, emphatic messages
  ###
  ### experimental cAIC minimization (completely experimental) # requires an input 'res' for .get_info_crits()
  if (FALSE && identical(processed$return_only,"cAICAPHLs")) {
    APHLs <- .get_info_crits(res)["cAIC"]
    if ( ! is.null(oldcAIC <- processed$port_env$objective)) {
      .update_port_fit_values(old_obj= - oldcAIC,new_obj= - APHLs[["cAIC"]], loopout_blob=loopout_blob,
                              models=models, processed=processed, control.HLfit=control.HLfit)
    } else if ( ! identical(control.HLfit$write_port_env,FALSE)) {
      processed$port_env$objective <- APHLs[["cAIC"]] # when return only cAIC => not fo residModel
    }
    res <- list(APHLs=APHLs)
    return(res)    ########################   R E T U R N
  }
  
  class(res) <- c("HLfit",class(res)) 
  # cleanup: for diagnostic, use
  # sort(sapply(ls(<object>$envir), function(x)
  # +             object.size(get(x, envir = <object>$envir))),decreasing=TRUE)
  lsv <- c("lsv",ls())
  rm(list=setdiff(lsv,"res")) ## empties the whole local envir except the return value
  return(res)
  
}
