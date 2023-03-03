.loop_while_TRUE <- function(processed, loopin_blob, loopout_blob, warningEnv, 
                             iter=0L, info_for_conv_rC=NULL, prev_lik=-Inf, conv_logL=NA,
                             conv.lambda=FALSE, conv.corr=FALSE, conv.phi=FALSE,
                             #
                             w.resid=loopout_blob$w.resid,
                             muetablob=loopout_blob$muetablob,
                             wranefblob=loopout_blob$wranefblob,
                             u_h=loopout_blob$u_h,
                             v_h=loopout_blob$v_h,
                             corr_est=loopout_blob$corr_est,
                             beta_eta=loopout_blob$beta_eta,
                             LMatrices=loopout_blob$LMatrices,
                             lambda_est=loopout_blob$lambda_est,
                             phi_est=loopout_blob$phi_est,
                             phiBLOB=loopout_blob$phiBLOB, # something strange: this line has been missing yet 
                                                          # neither the tests (&long &extras...) nor R CMD check say a problems
                             #
                             mu=muetablob$mu,
                             #
                             nrand=loopin_blob$nrand,
                             n_u_h=loopin_blob$n_u_h,
                             ZAL=loopin_blob$ZAL,
                             ranCoefs_blob=loopin_blob$ranCoefs_blob,
                             intervalInfo=loopin_blob$intervalInfo,
                             pforpv=loopin_blob$pforpv,
                             maxit.mean=loopin_blob$maxit.mean,
                             etaFix=loopin_blob$etaFix,
                             ranFix=loopin_blob$ranFix,
                             phi.Fix=loopin_blob$phi.Fix,
                             init.HLfit=loopin_blob$init.HLfit,
                             control.HLfit=loopin_blob$control.HLfit,
                             verbose=loopin_blob$verbose,
                             spaMM_tol=loopin_blob$spaMM_tol,
                             need_ranefPars_estim=loopin_blob$need_ranefPars_estim,
                             lam_fix_or_outer_or_NA=loopin_blob$lam_fix_or_outer_or_NA,
                             need_simple_lambda=loopin_blob$need_simple_lambda,
                             which_inner_ranCoefs=loopin_blob$which_inner_ranCoefs,
                             whichAPHLs=loopin_blob$whichAPHLs,
                             LMMbool=loopin_blob$LMMbool,
                             max.iter=loopin_blob$max.iter,
                             std_dev_res_needed_4_inner_estim=loopin_blob$std_dev_res_needed_4_inner_estim,
                             off=loopin_blob$off,
                             #
                             y=processed$y,
                             cum_n_u_h=processed$cum_n_u_h,
                             models=processed$models,
                             prior.weights=processed$prior.weights
) {
  
  while ( TRUE ) {
    ##the main loop with steps: new linear predictor, new leverages, new phi, new w.resid, new lambda, new fn(lambda)
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
                                wranefblob=wranefblob, u_h=u_h, v_h=v_h, w.resid=w.resid, 
                                H_w.resid=H_w.resid, # undefined if nrand=0L 
                                phi_est=phi_est, pforpv=pforpv, verbose=verbose,
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
    } else {
      PHIblob <- NULL
      conv.phi <- TRUE
    } ## there is a phi.Fix, already -> phi_est
    ###############################################################
    ######### Dispersion Estimates for lambda #####################
    ###############################################################
    
    if (need_ranefPars_estim) { ## lambda must be estimated 
      levlam_bound <- 1 - .spaMM.data$options$regul_lev_lambda
      if (warningEnv$leveLam1 <- any(whichleveLam1 <- (leverages$ranef > levlam_bound))) { 
        leverages$ranef[whichleveLam1] <- levlam_bound
        warningEnv$whichleveLam1 <- whichleveLam1 
      } 
      ################## ranefEstargs must contain arguments for makeCoveEst1 => its names are constrained
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
      calcRanefPars_blob <- .calcRanefPars(HLfit_corrPars=init.HLfit$corrPars,
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
        APHLs_ZX <- .calc_APHLs_from_ZX(auglinmodblob, which=whichAPHLs, processed) # with old lambda est
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
        corr_est <- list(rho=.getPar(calcRanefPars_blob$HLfit_corrPars,"rho")) ## not spaMM 3.0 : single inner-optimized rho
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
  loopout_blob$w.resid <- w.resid
  loopout_blob$prior.weights <- prior.weights
  loopout_blob$muetablob <- muetablob
  loopout_blob$wranefblob <- wranefblob
  loopout_blob$u_h <- u_h
  loopout_blob$v_h <- v_h
  loopout_blob$corr_est <- corr_est
  loopout_blob$beta_eta <- beta_eta
  loopout_blob$conv.phi <- conv.phi
  loopout_blob$LMatrices <- LMatrices
  loopout_blob$lambda_est <- lambda_est
  loopout_blob$phi_est <- phi_est

  loopout_blob$conv.phi <- conv.phi
  loopout_blob$conv.lambda <- conv.lambda
  loopout_blob$conv.corr <- conv.corr
  loopout_blob$iter <- iter
  loopout_blob$conv_logL <- conv_logL
  loopout_blob$APHLs <- APHLs
  if (need_ranefPars_estim) loopout_blob$calcRanefPars_blob <- calcRanefPars_blob
  if (std_dev_res_needed_4_inner_estim) loopout_blob$leverages <- leverages
  loopout_blob$auglinmodblob <- auglinmodblob  
  loopout_blob$PHIblob <- PHIblob
  
}


.add_unscaled_X.pv_fixef <- function(res, processed, beta_eta, etaFix, X.pv=processed$AUGI0_ZX$X.pv) {
  if ( ! is.null(attr(X.pv,"scaled:scale"))) {
    beta_eta <- .unscale(beta=beta_eta, X=X.pv)
    res$X.pv <- .unscale(X.pv) ## lvalue usefully not in an environment
  } else {
    res$X.pv <- X.pv
  }
  res$fixef <- .calc_full_fixef(processed, beta_eta, etaFix) # results with possible NA's
  res
}

.wrap_wrap_SEM <- function(processed, ZAL, beta_eta, off, corr_est, init.lambda, lam_fix_or_outer_or_NA, 
                           LMatrices, verbose, BinomialDen, phi_est) {
  ## specif probit
  processed$SEMargs$qr_X <- qr(processed$AUGI0_ZX$X.pv) 
  init.lambda <- .fill_init_lambda(init.lambda=init.lambda, processed=processed, ZAL=ZAL, cum_n_u_h=processed$cum_n_u_h)
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
  w.resid <- .calc_w_resid(muetablob$GLMweights, phi_est, obsInfo=processed$how$obsInfo) ## 'weinu', must be O(n) in all cases
  wranefblob <- processed$updateW_ranefS(u_h=u_h,v_h=v_h,lambda=lambda_est) 
  loopout_blob <- list(u_h=u_h, v_h=v_h, APHLs=APHLs, muetablob=muetablob, w.resid=w.resid, wranefblob=wranefblob,
                       iter=0L, conv.phi=FALSE, conv.lambda=FALSE, conv.corr=FALSE, conv_logL=NA, 
                       SEMblob=SEMblob, # not found in other loopout_blob's
                       corr_est=corr_est, prior.weights=processed$prior.weights)
  loopout_blob
  ##
}

.maxit.mean <- function(nrand, pforpv, etaFix, LMMbool, intervalInfo, inner_est_disp_pars, processed, models, phi.Fix) {
  if (nrand) { # (models[["eta"]]=="etaHGLM") {
    if (pforpv==0L && !is.null(etaFix$v_h)) {
      maxit.mean <- 0L ## used in test near the end...
    } else if ( LMMbool && is.null(intervalInfo) ) {
      maxit.mean <- 1L ## 1 is sufficient for LMM as Hessian does not vary with beta_eta  => quadratic function;
      # ./ E X C E P T for calculation of confidence intervals: at least two intervalSteps are required. Quite visible when 
      # dispersion params are outer-estimated (fitme) in which case there isno outer iteration to compensate a small maxit.mean  
    } else { ## even h maximization in *G*LMMs 
      if ( ! inner_est_disp_pars) { ## "hence no true outer iteration" (maybe a rough, outdated interpretation) 
        maxit.mean <- processed$iter_mean_dispFix 
      } else maxit.mean <- processed$iter_mean_dispVar # If phi.Fix and lam_fix_or_outer_or_NA, the only way to have 'outer' convergence is to have 'inner' convergence
    } 
  } else if (models[[1L]] %in% c("etaGLM")) {
    if ( ! .anyNULL(phi.Fix)) { ## 
      maxit.mean <- processed$iter_mean_dispFix 
    } else maxit.mean <- processed$iter_mean_dispVar # If phi.Fix and lam_fix_or_outer_or_NA, the only way to have 'outer' convergence is to have 'inner' convergence
  }
  maxit.mean
}

.add_loopout_vars <- function(loopout_blob, phi_est, inits_by_xLM, vec_nobs, eta, BinomialDen, processed, nrand, models, 
                              init.lambda, ZAL, cum_n_u_h, vec_n_u_h, n_u_h, ranCoefs_blob, etaFix, init.HLfit) {
  if ( ! is.null(vec_nobs)) {
    loopout_blob$PHIblob <- list(multiPHI=rep(list(list()), length(vec_nobs))) # consistently with .calcmultiPHIs() output 
  } else loopout_blob$PHIblob <- list()
  ## Finalize initial values for phi 
  loopout_blob$phi_est <- phi_est <- .denullify(phi_est, modifier=inits_by_xLM$phi_est, vec_nobs=vec_nobs)
  loopout_blob$muetablob <- muetablob <- .muetafn(eta=eta,BinomialDen=BinomialDen,processed=processed, phi_est=phi_est) 
  # with muetablob$mu being if Bin/Pois, O(n): facteur BinomialDen dans la transfo mu -> eta ## nonstandard mu des COUNTS
  loopout_blob$w.resid <- .calc_w_resid(muetablob$GLMweights, phi_est, obsInfo=processed$how$obsInfo) ## 'weinu', must be O(n) in all cases
  if (nrand) { # (models[["eta"]]=="etaHGLM") {
    ## Finalize initial values for lambda
    loopout_blob$lambda_est <- .HLfit_finalize_init_lambda(models, init.lambda, processed, ZAL, cum_n_u_h, 
                                                           vec_n_u_h, n_u_h, ranCoefs_blob) # mv __FIXME__ ->.calc_fam_corrected_guess() uses total nrand rather than nrand for submodels that contain the ranef
    loopout_blob$u_h <- processed$u_h_v_h_from_v_h(loopout_blob$v_h, lower.v_h=NULL, upper.v_h=NULL)
    loopout_blob$wranefblob <- processed$updateW_ranefS(u_h=loopout_blob$u_h,v_h=loopout_blob$v_h,lambda=loopout_blob$lambda_est) ## initialization !
  }
  # loopout_blob # envir, no return needed
}

