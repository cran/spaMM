HLfit_body_old <- function(processed, 
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
  prior.weights <- processed$prior.weights
  
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
  HLfit_corrPars <- init.HLfit$corrPars 
  if (! is.null(corr_est <- .get_cP_stuff(init.HLfit,"rho"))) corr_est <- list(rho=corr_est) 
  if (need_simple_lambda <- need_ranefPars_estim <- (nrand>0L)) {
    need_simple_lambda <- any(is.na(lam_fix_or_outer_or_NA) & ! ranCoefs_blob$is_set)
    need_ranefPars_estim <-  (need_simple_lambda || ! is.null(corr_est))
  } 
  #
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
  
  ### case where nothing to fit #############################################
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
  
  ### Initial values  for lambda, phi and beta from lam_fix_or_outer_or_NA, phi.Fix, or init.HLfit ##### 
  ## Initial estimate for phi 
  vec_nobs <- processed$vec_nobs
  phi_est <- phi.Fix
  phi_est <- .denullify(phi_est, modifier=processed$port_env$port_fit_values$phi_est, vec_nobs=vec_nobs)
  phi_est <- .denullify(phi_est, modifier=init.HLfit$phi, vec_nobs=vec_nobs)
  if (nrand) { 
    #
    if (identical(processed$return_only,"p_vAPHLs")) {
      whichAPHLs <- "p_v"
    } else if (identical(processed$return_only,"p_bvAPHLs")) {
      whichAPHLs <- "p_bv" ## return value may still include p_v if it is used to compute p_bv
    } else whichAPHLs <- c("p_v","p_bv")
    #
    ## Initial estimate for u_h, v_h 
    v_h <- processed$intervalInfo$init_v_h
    if (is.null(v_h)) v_h <- .initialize_v_h(processed, etaFix=etaFix, init.HLfit=init.HLfit) ## checks init.HLfit$v_h
    u_h <- processed$u_h_v_h_from_v_h(v_h, lower.v_h=NULL, upper.v_h=NULL)
    ## Initial estimate for lambda in 'compact" form
    init.lambda <- .calc_initial_init_lambda(lam_fix_or_outer_or_NA, nrand, processed, ranCoefs_blob, init.HLfit, ranFix)
  } else {
    u_h <- v_h <- numeric(0)
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
  intervalInfo <- processed$intervalInfo
  if (!is.null(intervalInfo)) {
    parmcol <- attr(intervalInfo$parm,"col")
    beta_eta[parmcol] <- intervalInfo$init ## already appropriately scaled if X.pv has been scaled
  }  
  ## Finalize initial values for phi 
  phi_est <- .denullify(phi_est, modifier=inits_by_xLM$phi_est, vec_nobs=vec_nobs)
  ## Finalize initial values for lambda
  if (nrand) { # (models[["eta"]]=="etaHGLM") {
    lambda_est <- .HLfit_finalize_init_lambda(models, init.lambda, processed, ZAL, cum_n_u_h, 
                                              vec_n_u_h, n_u_h, ranCoefs_blob) # mv __FIXME__ ->.calc_fam_corrected_guess() uses total nrand rather than nrand for submodels that contain the ranef
  }
  ###
  ## predictor from initial values. When there an etaFix, it is here in the 'off'set, and the dims of X.pv and beta_eta are here correspondingly reduced.
  if (nrand) { # (models[["eta"]]=="etaHGLM") {
    eta <- off + drop(processed$AUGI0_ZX$X.pv %*% beta_eta) + drop(ZAL %id*% v_h)
  } else  eta <- off +drop(processed$AUGI0_ZX$X.pv %*% beta_eta) ## no iteration hence no updating  ## FREQS
  ## conversion to mean of response variable (COUNTS for binomial)
  muetablob <- .muetafn(eta=eta,BinomialDen=BinomialDen,processed=processed, phi_est=phi_est) 
  mu <- muetablob$mu ## if Bin/Pois, O(n): facteur BinomialDen dans la transfo mu -> eta ## nonstandard mu des COUNTS
  w.resid <- .calc_w_resid(muetablob$GLMweights, phi_est, obsInfo=processed$how$obsInfo) ## 'weinu', must be O(n) in all cases
  conv.phi <- conv.lambda <- conv.corr <- FALSE
  std_dev_res_needed_4_inner_estim <- .anyNULL(phi.Fix) || need_simple_lambda 
  if (nrand) { # (models[["eta"]]=="etaHGLM") {
    if (! is.null(intervalInfo)) intervalInfo$ranFix <- ranFix
    wranefblob <- processed$updateW_ranefS(u_h=u_h,v_h=v_h,lambda=lambda_est) ## initialization !
    if (pforpv==0L && !is.null(etaFix$v_h)) {
      maxit.mean <- 0L ## used in test near the end...
    } else if ( LMMbool && is.null(intervalInfo) ) {
      maxit.mean <- 1L ## 1 is sufficient for LMM as Hessian does not vary with beta_eta  => quadratic function;
      # ./ E X C E P T for calculation of confidence intervals: at least two intervalSteps are required. Quite visible when 
      # dispersion params are outer-estimated (fitme) in which case there isno outer iteration to compensate a small maxit.mean  
    } else { ## even h maximization in *G*LMMs 
      # for inner ranCoef estim [eg, HLfit3 example with family=Gamma(log)] no leverages needed => std_dev_res_needed_4_inner_estim is FALSE 
      # yet if there are inner_ranCoefs and only inner estim, not outer, is used, iter_mean_dispVar seems appropriate. 
      # => second condition on which_inner_ranCoefs: [not a strict check of no outer estime but probably equivalent in practice]
      if ( ! (std_dev_res_needed_4_inner_estim || length(which_inner_ranCoefs))) { 
        maxit.mean <- processed$iter_mean_dispFix 
      } else maxit.mean <- processed$iter_mean_dispVar # If phi.Fix and lam_fix_or_outer_or_NA, the only way to have 'outer' convergence is to have 'inner' convergence
    } 
  } else if (models[[1L]] %in% c("etaGLM")) {
    if (! is.null(intervalInfo)) {
      intervalInfo$parmcol_X <- parmcol 
    }
    if ( ! .anyNULL(phi.Fix)) { ## 
      maxit.mean <- processed$iter_mean_dispFix 
    } else maxit.mean <- processed$iter_mean_dispVar # If phi.Fix and lam_fix_or_outer_or_NA, the only way to have 'outer' convergence is to have 'inner' convergence
  }
  prev_lik <- -Inf
  conv_logL <- NA
  ########################################
  ######### Main loop ####################
  ########################################
  iter <- 0L
  if ( ! is.null(vec_nobs)) {
    PHIblob <- list(multiPHI=rep(list(list()), length(vec_nobs))) # consistently with .calcmultiPHIs() output 
  } else PHIblob <- list()
  info_for_conv_rC <- NULL
  if (HL[1]=="SEM") { ## specif probit
    processed$SEMargs$qr_X <- qr(processed$AUGI0_ZX$X.pv) 
    init.lambda <- .fill_init_lambda(init.lambda=init.lambda, processed=processed, ZAL=ZAL, cum_n_u_h=cum_n_u_h)
    locarglist <- list(processed=processed, ZAL=ZAL, beta_eta=beta_eta,
                       off=off, corr_est=corr_est, init.lambda=init.lambda,
                       lambda.Fix=lam_fix_or_outer_or_NA, LMatrices=LMatrices, verbose=verbose)
    SEMblob <- .probitgemWrap("SEMwrap",arglist=locarglist, pack="probitgem") # eval(as.call(c(quote(SEMwrap),logarglist)))
    beta_eta <- SEMblob$beta_eta
    corr_est["rho"] <- .get_rho_from_SEM_output(SEMblob, lam_fix_or_outer_or_NA)
    if (is.null(SEMblob$glm_lambda)) {
      lambda_est <- lam_fix_or_outer_or_NA
    } else {
      lambda_est <- predict(SEMblob$glm_lambda,type="response")
    }
    u_h <- v_h <- SEMblob$v_h
    logLapp <- SEMblob$logLapp
    attr(logLapp,"seInt") <- SEMblob$seInt ## may be NULL
    APHLs <- list(logLapp=logLapp) ## keeps attributes
    ## for summary.HLfit -> beta_cov (_info ?)
    eta <- off + drop(processed$AUGI0_ZX$X.pv %*% beta_eta) + drop(ZAL %id*% v_h)
    muetablob <- .muetafn(eta=eta,BinomialDen=BinomialDen,processed=processed, phi_est=phi_est) 
    APHLs$clik <- .calc_clik(muetablob$mu,phi_est,processed) ## useful for .get_info_crits(); phi_est must be trivial an not updated below.
    w.resid <- .calc_w_resid(muetablob$GLMweights,phi_est, obsInfo=processed$how$obsInfo) ## 'weinu', must be O(n) in all cases
    wranefblob <- processed$updateW_ranefS(u_h=u_h,v_h=v_h,lambda=lambda_est) 
    ##
  } else while ( TRUE ) { ##the main loop with steps: new linear predictor, new leverages, new phi, new w.resid, new lambda, new fn(lambda)
    if (nrand) { # (models[["eta"]]=="etaHGLM") {
      H_w.resid <- .calc_H_w.resid(w.resid, muetablob=muetablob, processed=processed) # for LLF w.resid is not generally defined.
      if (! is.null(intervalInfo)) { # There is a distinct interval procedure in .calc_etaGLMblob()
        intervalInfo$corr_est <- corr_est ## currently needs to be added ex tempo... apparently only for warn_intervalStep()
        intervalInfo$beta_eta <- beta_eta ## ex tempo: the current estimate of the CI bound
        # to compute a likelihood for theta:=(beta, lambda, phi_pars) The predictions of any ranef must be the ones under theta.
        # This means that if these predictions change with beta the lik does not really be computed 
        # unless, beta,V have reached equilibrium (the present end of the loop).
        # However, in case of mixed model for phi, we would need to *predict* the new phi values as function of given phi parameters
        # and of new deviance residuals implied by the (beta, lambda) and new (beta,V) 
        # (for a phi GLM, no such problem: the phi pars determine the phi values irrespective of the deviance residuals)
        # Hence we can warn about suspiciously high likelihood when:
        intervalInfo$phi_pred_OK <- (intervalInfo$no_phi_pred || conv.phi) 
        #where no_phi_pred is a perhaps too strict condition and || conv.phi possible a too loose one.
      } ## else intervalInfo remains NULL
    } 
    auglinmodblob <- .wrap_IRLS(nrand=nrand, intervalInfo=intervalInfo, processed=processed, beta_eta=beta_eta,
                                ZAL=ZAL, y=y, lambda_est=lambda_est, muetablob=muetablob, maxit.mean=maxit.mean, etaFix=etaFix, 
                                wranefblob=wranefblob, u_h=u_h, v_h=v_h, w.resid=w.resid, H_w.resid=H_w.resid, phi_est=phi_est, pforpv=pforpv, verbose=verbose,
                                ad_hoc_corrPars= {  # evaluated only for spprec.
                                  adj_rho <- corr_est$rho
                                  if (is.null(adj_rho)) adj_rho <- .getPar(ranFix,"rho") 
                                  if (is.null(adj_rho)) adj_rho <- init.HLfit$rho
                                  list(rho=adj_rho)
                                })
    if ( ! is.null(auglinmodblob)) { # may be NULL if formula was ~ 0 (or ~ offset, presumably)
      if ( is.null(processed$X_off_fn)) beta_eta <- auglinmodblob$beta_eta
      muetablob <- auglinmodblob$muetablob 
      mu <- muetablob$mu ## necess dans test COMPoisson HLfit...
      w.resid <- auglinmodblob$w.resid 
      if (nrand) {
        v_h <- auglinmodblob$v_h
        u_h <- auglinmodblob$u_h
        wranefblob <- auglinmodblob$wranefblob # updated only if !LMM
      }
    } else auglinmodblob <- list(muetablob=muetablob, w.resid=w.resid)
    
    # 'mu' must be a vector not a matrix. This was checked herefrom version 3.13 (at least) to version 4.1.66
    # Then I moved the check after the fit.

    ########## LEVERAGES
    if (std_dev_res_needed_4_inner_estim) {
      leverages <- .calc_std_leverages(models, need_ranefPars_estim=need_ranefPars_estim, phi.Fix=phi.Fix, auglinmodblob=auglinmodblob, 
                                       n_u_h=n_u_h, nobs=length(y), processed=processed, 
                                       # w.resid=w.resid, u_h=u_h, 
                                       need_simple_lambda=need_simple_lambda, 
                                       #muetablob=muetablob, wranefblob=wranefblob, 
                                       ZAL=ZAL, lambda_est=lambda_est, cum_n_u_h=cum_n_u_h, phi_est=phi_est)
    }
    ######### Dispersion Estimates for phi #####################
    if (.anyNULL(phi.Fix)) { ## if phi is estimated (vs phi.Fix set to 1 for Poisson, Binomial); implies std_dev_res_needed_4_inner_estim
      leverages$resid[leverages$resid>1-1e-8] <- 1-1e-8
      PHIblob <- .calcPHIs(processed=processed, 
                           y=y,mu=mu, wt= eval(prior.weights), 
                           lev_phi=leverages$resid, phimodels=models[["phi"]], verbose=verbose, 
                           control.HLfit=control.HLfit, phi.Fix=phi.Fix,
                           iter=iter, prev_PHIblob=PHIblob)

      .addPhiGLMwarning(PHIblob, phimodels=models[["phi"]], warningEnv)
      #if ( is.null(processed$data$.phi)) browser() # we're in a mean-response fit
      
      # phi_est may be (a single value|a long vector| a list of such values, for mv).
      conv.phi <- .eval_conv.phi(phi_est, next_phi_est=PHIblob$next_phi_est, # value of *phi* (not phi_i:= phi/prior.weights as pw are used in GLMweights, not here)  
                                 spaMM_tol=spaMM_tol)
      # :           # may substract scalar (initial value) to vector (predictions from phi model) 
    } else {conv.phi <- TRUE} ## there is a phi.Fix, already -> phi_est
    ###############################################################
    ######### Dispersion Estimates for lambda #####################
    ###############################################################
    
    if (need_ranefPars_estim) { ## lambda must be estimated 
      levlam_bound <- 1 - .spaMM.data$options$regul_lev_lambda
      if (warningEnv$leveLam1 <- any(whichleveLam1 <- (leverages$ranef > levlam_bound))) { 
        leverages$ranef[whichleveLam1] <- levlam_bound
        warningEnv$whichleveLam1 <- whichleveLam1 
      } 
      ################## ranefEstargs mustcontain arguments for makeCoveEst1 => its names are constrained
      ####
      ranefEstargs <- list(u_h=u_h,ZAlist=processed$ZAlist,cum_n_u_h=cum_n_u_h,
                           prev_LMatrices=LMatrices,
                           processed=processed,
                           init_ranCoefs=init.HLfit$ranCoefs,
                           H_w.resid=H_w.resid,
                           w.resid=w.resid,
                           prev_lambda_est=lambda_est)
      if (any(ranCoefs_blob$isRandomSlope)) { ## if random-slope model
        ranefEstargs <- c(ranefEstargs,list(phi_est=phi_est,
                                            as_matrix=( ! inherits(ZAL,"Matrix")),v_h=v_h))
        ## MakeCovEst defines à local ZAL and the eta,mu, w.resid must generally be recomputed locally for this ZAL
        # Thi is a list of arguments for .makeCovEst1 -> objfn -> .solve_IRLS_as_ZX  
        ranefEstargs$MakeCovEst_pars_not_ZAL_or_lambda <- list(
          muetablob=muetablob,
          maxit.mean=maxit.mean, etaFix=etaFix,
          ## supplement for LevenbergM
          beta_eta=beta_eta,
          ## supplement for ! GLMM
          u_h=u_h, v_h=v_h, phi_est=phi_est,
          for_intervals=intervalInfo, # inner-estimation of ranCoefs is incompatible with the confint hack.  
          #   Cf early detection in .confint_LRT(). If we bypass this early catch, setting for_intervals=NULL here
          # yields smaller innacuracies than setting it to intervalInfo, which is totally off. 
          # In a a non-confint fit It needs to be set to (intervalInfo=NULL) so that the solve_IRLS function will know it's non-confint. 
          verbose=c(TRACE=FALSE,trace=FALSE),
          ##
          processed=processed)
      }
      calcRanefPars_blob <- .calcRanefPars(HLfit_corrPars=HLfit_corrPars,
                                           lev_lambda=leverages$ranef,
                                           ranefEstargs=ranefEstargs,
                                           ranCoefs_blob=ranCoefs_blob,
                                           lam_fix_or_outer_or_NA=lam_fix_or_outer_or_NA,
                                           rand.families=processed$rand.families,
                                           psi_M=processed$psi_M,
                                           verbose=verbose,
                                           iter=iter,
                                           control=processed[["control.glm"]],
                                           maxLambda=processed$maxLambda,
                                           warningEnv)
      ####
      HLfit_corrPars <- calcRanefPars_blob$HLfit_corrPars
      ########################################
      # => variation of log(u^2/lambda) = simplified likRanU convergence  (from 1.9.31) (see below use of verbose["trace"] to print diagnostics)
      # conv lambda code depending on gaussian_u_ranges had been ineffective for a long time and commented-out in v. 3.6.38 
      # then fully removed in v.3.6.39. Also ***conv.lambda no longer depends on iter>1L***.
      conv.lambda <- .eval_conv.lambda(next_lambda_est=calcRanefPars_blob$next_lambda_est, 
                                       phi.Fix=phi.Fix, phi_est=phi_est, spaMM_tol=spaMM_tol, lambda_est=lambda_est)
      lambda_est <- calcRanefPars_blob$next_lambda_est
      ## convergence of corr estimate
      conv.corr <- .eval_conv.corr(next_info_for_conv_rC=calcRanefPars_blob$HLfit_corrPars$info_for_conv_rC, ## vector 
                                   iter=iter, info_for_conv_rC=info_for_conv_rC, spaMM_tol=spaMM_tol)
    } else { conv.lambda <- conv.corr <- TRUE } ## end if need_ranefPars_estim else...
    #
    iter <- iter+1L ## here first from 0 to 1
    ###### convergence: 
    if ( (conv.phi && conv.lambda && conv.corr) ||
         iter==max.iter ## did not converge...
         ) {
      ## the final APHLs should always be with updated lambda est (makes a difference in "hands-on calculations" in the gentle intro) ...
      APHLs <- .update_APHLs(need_ranefPars_estim, muetablob, phi_est, lambda_est, u_h, v_h, ZAL, processed, whichAPHLs, 
                             corr_est, ranFix, init.HLfit, auglinmodblob, nrand)
      if ( ! need_ranefPars_estim || iter== max.iter) {
        break
      } else {
        ## => perform a stricter check of possible convergence of logL 
        ## Without this check , taking final APHLs from .calc_APHLs_from_params() rather than from .calc_APHLs_from_ZX() changed the test-Rasch final lik. 
        ##    Also I would see greater inaccuracies in the "independent fits" mv tests.
        APHLs_ZX <- .calc_APHLs_from_ZX(auglinmodblob,which=whichAPHLs,processed) # with old lambda est
        if (max(abs(APHLs[[processed$objective]]-APHLs_ZX[[processed$objective]]))<
            .spaMM.data$options$spaMM_tol$logL_tol) break
        # else there is too much variation in processed$objective: the loop continues.
      }
    } 
    ##
    ## Perform various updates and continue (but one possible break after several updates)
    if (.anyNULL(phi.Fix)) phi_est <- PHIblob$next_phi_est # value of *phi* (not phi_i:= phi/prior.weights as pw are used in GLMweights, not here)
    if (nrand) { # (models[["eta"]]=="etaHGLM") {
      if ( ! is.null(corr_est) ) {
        corr_est <- list(rho=.getPar(HLfit_corrPars,"rho")) ## not spaMM 3.0 : single inner-optimized rho
        #LMatrix est constant!= decomp$u
      }
      if (need_ranefPars_estim) { ## next_lambda_est/ next ranCoefs/ next rho adjacency available:
        if (any(ranCoefs_blob$isRandomSlope)) {
          LMatrices <- calcRanefPars_blob$next_LMatrices ## keep L factor of corr mats for all models 
          if (processed$is_spprec) {
            for (rt in which_inner_ranCoefs) {
              processed$AUGI0_ZX$envir$precisionFactorList[[rt]] <- 
                .precisionFactorize(latentL_blob=attr(LMatrices[[rt]],"latentL_blob"),
                                    rt=rt, longsize=ncol(LMatrices[[rt]]), processed=processed,
                                    cov_info_mat=processed$corr_info$cov_info_mats[[rt]])
            }
          }
          ZAL <- .compute_ZAL(XMatrix=LMatrices, ZAlist=processed$ZAlist,as_matrix=( ! inherits(ZAL,"Matrix")),
                              bind.= ! processed$is_spprec) 
          if ( ! LMMbool ) {
            ## ZAL is modified hence wranefblob must be modified (below) but also eta-> mu->GLMweights
            ## .makeCovEst1 may have reestimated beta but we do not take this into account nor any resulting change in the 'blobs'
            eta <- off + drop(processed$AUGI0_ZX$X.pv %*% beta_eta) + drop(ZAL %id*% v_h) 
            muetablob <- .muetafn(eta=eta,BinomialDen=processed$BinomialDen,processed=processed, phi_est=phi_est) 
          }
        }
        # UPDATE:
        ##      in particular, also in the adjacency case if rho was updated but not the lambda param.
        wranefblob <- processed$updateW_ranefS(u_h=u_h,v_h=v_h,lambda=lambda_est) ## bc lambda was modified
        info_for_conv_rC <- calcRanefPars_blob$HLfit_corrPars$info_for_conv_rC ## vector
      } 
    }
    if ( .anyNULL(phi.Fix) || ( nrand && any(ranCoefs_blob$isRandomSlope) && ! LMMbool) ) { ## phi or (ZAL -> mu) modified
      w.resid <- .calc_w_resid(muetablob$GLMweights,phi_est, obsInfo=processed$how$obsInfo) ## bc phi was updated. 'weinu', must be O(n) in all cases 
    }
    ## conv_logL either used to break the loop, Xor required only in last two iters for diagnostics 
    if (processed$break_conv_logL || verbose["trace"] ) {
      next_lik <- .calc_APHLs_from_ZX(auglinmodblob=auglinmodblob, which="p_v", processed=processed)$p_v
      # this does not depend on the latest ranPars estimates, since sXaug was not updated afterwards... 
      conv_logL <- abs(next_lik - prev_lik)/(0.1 + abs(next_lik)) < 1e-8 # ~ glm.fit convergence
      if (processed$break_conv_logL && conv_logL) break 
      prev_lik <- next_lik
    } else conv_logL <- NA
    ##
    if (verbose["trace"]) {
      print(paste("iteration ",iter,
                  #"; convergence criteria for phi, lambda, corr pars , conv_lambda_vs_u, conv_rel_lambda: ",
                  "; convergence criteria for phi, lambda, corr pars: ",
                  paste0(c( conv.phi , conv.lambda, conv.corr), #,conv_lambda_vs_u, conv_rel_lambda),
                         collapse = " "),
                  "; logL & conv_logL:", next_lik, conv_logL))
    } 
    ##### end convergence block
  } ## end main loop while ( TRUE )
  ########################################
  ######### END main loop ################
  ########################################
  if (verbose["trace"]) {
    if (iter==max.iter) {
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
      res <- list(APHLs=APHLs) 
      if ( ! is.null(oldp_v <- processed$port_env$objective)) {
        .update_port_fit_values(old_obj=oldp_v,new_obj=APHLs$p_v, 
                                port_fit_values=list(fixef=beta_eta,v_h=v_h,phi_est=phi_est), 
                                models=models, processed=processed, control.HLfit=control.HLfit,
                                lambda_est=lambda_est, 
                                PHIblob=PHIblob)
      } else if ( ! identical(control.HLfit$write_port_env,FALSE)) {
        processed$port_env$objective <- APHLs$p_v # when return only p_v => not for residModel
      }
    }
    return(res)   ########################   R E T U R N
  } 
  ## ELSE continue: make sure p_bv is included
  if (HL[1] != "SEM") {
    if (nrand) { # (models[["eta"]]=="etaHGLM") {
      ## cf notes 19/08/2016 pour calcul APHLs et IC's for phiHGLM 
    } else  APHLs <- .calc_APHLs_XLM(processed, w.resid, clik=APHLs$clik) ## GLM|LLM. 'w.resid' used only for 'REML' which is a new concept for a fixed-effect LLM 
  }
  ######################### potential R E T U R N here: with p_bv
  if ( identical(processed$return_only,"p_bvAPHLs")) {
    if ( ! is.null(oldp_bv <- processed$port_env$objective)) {
      .update_port_fit_values(old_obj=oldp_bv,new_obj=APHLs$p_bv, 
                              port_fit_values=list(fixef=beta_eta,v_h=v_h,phi_est=phi_est), 
                              models=models, processed=processed, control.HLfit=control.HLfit,
                              lambda_est=lambda_est, 
                              PHIblob=PHIblob)
    } else if ( ! identical(control.HLfit$write_port_env,FALSE)) {
      processed$port_env$objective <- APHLs$p_bv # when return only p_bv => not for residModel
    }
    res <- list(APHLs=APHLs)
    return(res)    ########################   R E T U R N
  } 
  # beta_cov code removed from here in v1.9.24
  ######################
  ######################
  ######################
  ######################################
  ## BUILD full RETURN VALUE
  ######################################
  #
  ###################
  ## LIKELIHOODS
  ###################
  res <- list(APHLs=APHLs)
  ###################
  ## FIXEF, ETA, ... 
  ###################
  res$X.pv <- processed$AUGI0_ZX$X.pv
  if ( ! is.null(attr(res$X.pv,"scaled:scale"))) {
    beta_eta <- .unscale(beta=beta_eta, X=res$X.pv)
    res$X.pv <- .unscale(res$X.pv) ## usefully not in an environment
  } 
  res$fixef <- .calc_full_fixef(processed, beta_eta, etaFix) # results with possible NA's
  if (identical(processed$return_only,"confint_bound")) {
    return(res)    ########################   R E T U R N fixef + APHLs
  }
  # sane_eta with names
  res$eta <- .format_eta(eta=muetablob$sane_eta, data=processed$data) ## convenient for defining starting values... and also sometimes used by predict() # now kept in muetablob
  ####
  res$muetablob <- muetablob # directly as a lot of elements may be needed:
    # $dmudeta, eg for dvdlogphiMat in.get_logdispObject() 
    # $p0=muetablob$p0, for simulate 
    # Md3logcLdeta3, Md2logcLdeta2 for hatvalues() -> .hatvals2std_lev() -> . -> ..calc_dlW_deta()
    # $mv for and hatvalues.HLfit(). (and possible for mv simulate)
  ###################
  ## DATA
  ###################
  res$data <- processed$data
  if (is.null(processed$family) || processed$family$family=="binomial") { # null for mv case
    res$BinomialDen <- BinomialDen # we could put it in all cases...
  }
  res$y <- y ## counts for Pois/bin
  res$prior.weights <- prior.weights # with attrs, and possibly a call
  res$prior.weights[] <- eval(prior.weights) ## see Gamma()$simulate # eval the call but keep the attributes.
  ###################
  ## MODEL info
  ###################
  res$family <- processed$family
  res$families <- processed$families
  res$ranFix <- ranFix ## currently as a uniform template consistent with projected changes ; except that lamFix, phiFix info is now in lambda.object, etc
  .canonize_disp_envs(fitobject=res)
  ## If an outer optimizer has been called,
  #  "fix" and "outer" parameters are given these types by .get_refit_args() after the optimization call, 
  # then HLfit called again and we reach this point.
  #  This means ranFix gets its type from there *if* properly retained by .canonizeRanPars() 
  #  Then we add inner-optimized parameters, with "var" type added by .get_CorrEst_and_RanFix()
  if ( ! is.null(corr_est) && ! is.null(init.HLfit$corrPars)) corr_est <- list(corrPars=relist(corr_est$rho,init.HLfit$corrPars)) ## not yet spaMM 3.0
  # Canonical, and inherits all info about outer-optimized corrPars through HLfit's ranFix argument:
  res$CorrEst_and_RanFix <- .get_CorrEst_and_RanFix(ranFix, corr_est) # corr_est parameters are inner-estimated and of type "var"
  if ( ! is.null(res$CorrEst_and_RanFix$corrPars)) {
    res$corrPars <- structure(res$CorrEst_and_RanFix$corrPars, # ## subset of the above: F I X M E (?) redundancy but convenient when examining fits
                              type=attr(res$CorrEst_and_RanFix,"type")$corrPars,
                              message='Use get_ranPars(.,which="corrPars") to extract "corrPars" cleanly from fit object.')
  }  
  #
  ##### LAMBDA and other RANEF PARS
  if (need_ranefPars_estim) { # (FALSE if only outer_ranCoefs)
    process_resglm_blob <- .bloc_lambda(models=models, 
                                        processed=processed, lam_fix_or_outer_or_NA=lam_fix_or_outer_or_NA, 
                                        cum_n_u_h=cum_n_u_h, next_LMatrices=LMatrices,
                                        HL=HL,
                                        SEMblob=SEMblob, # only for SEM, versus next ones for non-SEM:
                                        calcRanefPars_blob=calcRanefPars_blob,
                                        lev_lambda=leverages$ranef)
  } else {
    process_resglm_blob <- list(lambda_pred_list=as.list(rep(NA,nrand)))
  }
  res$dfs <- .calc_dfs(need_ranefPars_estim, process_resglm_blob, init.lambda, 
                   ranCoefs_blob, LMatrices, processed, pforpv, CorrEst_and_RanFix=res$CorrEst_and_RanFix)
  res$models <- models # structure(models, LMMbool=LMMbool) # $models already has this attribute !
  res$main_terms_info <- processed$main_terms_info ## used by predict
  res$predictor <- processed$predictor ##  all post fitting functions expect PROCESSED predictor
  res$vec_nobs <- processed$vec_nobs ## non-null for fitmv
  #
  res$REMLformula <- processed$REMLformula  # only for .REMLmess()... but it's still a simple way to pass the info. Perhaps put it elsewhere in res?
  ###################
  ## OBJECTIVE and ALGORITHMs
  ###################
  res$HL <- HL ## info on fitting objective
  res$how <- list(spaMM.version=packageVersion("spaMM"),
                  MME_method=.get_MME_method(auglinmodblob, HL),
                  switches=c(augZXy_cond=processed$augZXy_cond,
                             use_spprec_QR=.spaMM.data$options$use_spprec_QR),
                  obsInfo=processed$how$obsInfo
  )
  # res$MME_method <- structure(res$how$MME_method,
  #                             message="Please use how(<fit object>)[['MME_method']] to extract this information cleanly.")                 
  res$spaMM.version <- structure(res$how$spaMM.version, ## this is NOT a string and comparison with a string is suitably def'ed (as detailed in ?package_version)
                                 message="Please use how(<fit object>)[['spaMM.version']] to extract this information cleanly.")                 
  if (HL[1]=="SEM") {
    res$SEM_info <- SEMblob$SEM_info ## info
  } # else  
  ###################
  ## FITTED VALUES
  ###################
  res$fv <- .mu_U2fv(mu, # For truncated fams, input mu of latent untruncated variable, and muetablob is source of additional info; output fv is different but mu is kept as an attribute
                     # For binomial: input mu is count, output fv is proba
                     BinomialDen=BinomialDen, muetablob=muetablob, processed=processed)
  ###################
  ## LEVERAGES and REML (ie either phi OR lambda was estimated)
  ###################
  if (HL[1]!="SEM") { ## both lev_phi and deviance_residual missing otherwise
    if (std_dev_res_needed_4_inner_estim) { ## ll model leverages are computed and it makes sense to consider the residuals
      res$lev_phi <- leverages$resid
      dev_res_blob <- .std_dev_resids(res, phi_est=phi_est, lev_phi=leverages$resid)
      if (inherits(mu,"Matrix")) {
        warning("inefficiency detected. Please contact the package maintainer.", immediate. = TRUE) # it's inefficient if true in the loop...
        mu <- drop(mu) ## Old comment: "pb calcul deviance_residual" which is why I moved the test from the main loop to here.
      }
      res$std_dev_res <- sign(y-mu) * dev_res_blob$std_dev_res
      dev_res <- dev_res_blob$dev_res # needed below
      
      if (need_simple_lambda) res$lev_lambda <- leverages$ranef # __F I X M E__ remove the local condition ?
      
      # res$diagnostics$m_grad_obj <- auglinmodblob$m_grad_obj # typically NULL for LMM
      if (nrand && is.null(warningEnv$leveLam1)) { # .calcRanefPars was not called
        # Cf singfitF test in test-rank: outer estim without refit => reaches this point
        #But we want to diagnose possible non-identifiability of lambda => add info to warningEnv:
        .diagnose_lev_lambda(leverages, warningEnv, nrand, cum_n_u_h)
      }
    }
  }  
  res$distinctX.Re <- processed$X.Re ## NULL if not distinct from X.pv
  ###################
  ## ALL other LAMBDA returns
  ###################
  res$rand.families <- processed$rand.families 
  ##
  res$QRmethod <- processed$QRmethod
  #
  ## $w.ranef and $w.resid not doc'ed, as there is no mention of the augmented model in the doc.
  # if (is.list(w.resid)) { ## truncated 'family', or 'families' with some truncated one(s), but not all mv cases
  #   res$w.resid <- w.resid$w_resid ## useful for .get_info_crits() and get_LSmatrix
  # } else res$w.resid <- w.resid ## useful for .get_info_crits() and get_LSmatrix()
  # # res$w.resid is always a vector
  # if (is.null(attr(res$w.resid,"unique"))) attr(res$w.resid,"unique") <- length(unique(res$w.resid))==1L # is.null() => for mv    
  # #
  if (nrand) { # (models[["eta"]]=="etaHGLM") {
    res <- .add_ranef_returns(res, processed=processed, wranefblob=wranefblob, lambda_est=lambda_est, 
                              process_resglm_blob=process_resglm_blob, LMatrices=LMatrices, init.lambda=init.lambda, 
                              v_h=v_h, u_h=u_h, ranCoefs_blob=ranCoefs_blob, moreargs=attr(fixed,"moreargs"))
  } ## else various res$ elements are NULL
  ###################
  ## ALL other PHI returns
  ###################
  res <- .add_phi_returns(res=res, processed=processed, PHIblob=PHIblob, phi.Fix=phi.Fix, dev_res=dev_res, phi_est=phi_est)
  ################### the magic environment
  res$envir <- .add_fitobject_envir(nrand=nrand, HL=HL, w.resid=w.resid, auglinmodblob=auglinmodblob, processed=processed)
  ###################
  ## WARNINGS
  ###################
  .hack_options_error(message=NULL)
  ## translation of warnings in user-more friendly form 
  res$warnings <- .post_process_warningEnv(warningEnv=warningEnv, processed=processed, maxit.mean=maxit.mean,  pforpv=pforpv, 
                                           innerj=auglinmodblob$innerj, nonSPD=.BLOB(res$envir$sXaug)$nonSPD, conv_logL=conv_logL, 
                                           iter=iter, conv.lambda=conv.lambda, conv.phi=conv.phi, conv.corr=conv.corr)
  .verbose_warnings(verbose, res$warnings) # immediate, emphatic messages
  ###
  ### experimental cAIC minimization (completely experimental)
  if (FALSE && identical(processed$return_only,"cAICAPHLs")) {
    APHLs <- .get_info_crits(res)["cAIC"]
    if ( ! is.null(oldcAIC <- processed$port_env$objective)) {
      .update_port_fit_values(old_obj= - oldcAIC,new_obj= - APHLs[["cAIC"]], 
                              port_fit_values=list(fixef=beta_eta, # '___F I X M E___' is the unscaled one, should be the scaled one... fixed in 'loopout' version
                                                   v_h=v_h,phi_est=phi_est), 
                              models=models, processed=processed, control.HLfit=control.HLfit,
                              lambda_est=lambda_est, 
                              PHIblob=PHIblob)
    } else if ( ! identical(control.HLfit$write_port_env,FALSE)) {
      processed$port_env$objective <- APHLs[["cAIC"]] # when return only cAIC => not fo residModel
    }
    res <- list(APHLs=APHLs)
    return(res)    ########################   R E T U R N
  }
  ###################
  ## SUMMARY, RETURN
  ###################
  class(res) <- c("HLfit",class(res)) 
  # cleanup: for diagnostic, use
  # sort(sapply(ls(<object>$envir), function(x)
  # +             object.size(get(x, envir = <object>$envir))),decreasing=TRUE)
  lsv <- c("lsv",ls())
  rm(list=setdiff(lsv,"res")) ## empties the whole local envir except the return value
  return(res)
}
